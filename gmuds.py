import os
import sys
import pickle
import datetime
import copy
from itertools import groupby
from operator import itemgetter
import glob

import numpy as np
import pandas as pd
import pyproj
import xarray as xr
from affine import Affine
import rasterio as rio

_PATH = os.path.dirname(os.path.realpath(__file__))
print('Source Path: {}'.format(_PATH))

sys.path.append(_PATH)

DAT_PATH = os.path.join(_PATH, 'dat', 'GMU_DS_product')
print('Data Path: {}'.format(DAT_PATH))

AMB_PATH = os.path.join(DAT_PATH, 'Ambient')
SWR_PATH = os.path.join(DAT_PATH, 'SW')
WND_PATH = os.path.join(DAT_PATH, "UVWind")
PRC_PATH = os.path.join(DAT_PATH, 'DRain_hyb')

TRISHULI_PATH = os.path.join(_PATH, 'dat', 'gmu_trishuli')

# class
class _Grid(object):

    def __init__(self):
        self.p = pyproj.Proj(self._p_str)
        self.gt = (self._x0, self._res, 0, self._y0, 0, -self._res)
        self.fwd = Affine.from_gdal(*self.gt)
        self.rev = ~self.fwd
        # to get coords
        self.x_size = self._res * self._shape[1]
        self.y_size = self._res * self._shape[0]

    @property
    def bounds(self):
        return (self._x0, self._x0 + self.x_size,
                self._y0, self._y0 - self.y_size)

    @property
    def coords(self):
        return [
            (self._x0, self._y0),
            (self._x0 + self.x_size, self._y0),
            (self._x0 + self.x_size, self._y0 - self.y_size),
            (self._x0, self._y0 - self.y_size),
            (self._x0, self._y0)
        ]

    def get_max_bounds(self, grid):
        x_, y_ = [], []
        for (x, y) in grid.coords:
            x_c, y_c = grid.transform(self, x, y)
            x_.append(x_c)
            y_.append(y_c)
        return min(x_), max(x_), min(y_), max(y_)

    def transform(self, grid, x, y):
        return pyproj.transform(self.p, grid.p, x, y)

    def itransform(self, grid, x, y):
        return pyproj.transform(grid.p, self.p, x, y)

    def forward(self, col, row, grid=None):
        """
        Get x, y map coordinates from col, row in *current* grid.

        :param col: column index
        :param row: row index
        :param grid: grid that contains coordinates
        :return: (x, y) - grid coordinates
        """
        x, y = self.fwd * (col, row)
        if grid is not None:
            x, y = self.transform(grid, x, y)
        return x, y

    def inverse(self, x, y, grid=None):
        """
        Get col, row indices from *any* map grid coordinates.
        :param x: scan position coordinate
        :param y: raster position coordinate
        :param grid: grid that contains coordinates
        :return: (col, row) - array indices
        """
        if grid is not None:
            x, y = grid.transform(self, x, y)
        col, row = self.rev * (x, y)
        return np.floor(col).astype('int'), np.ceil(row).astype('int')  #


class GmudsGrid(_Grid):
    _p_str = """+proj=lcc +lat_1=30 +lat_2=62 +lat_0=0 +lon_0=105 +x_0=0 +y_0=0 
                +ellps=WGS84 +datum=WGS84 +units=m +no_defs"""
    _res = 1000
    _x0 = -3000000
    _y0 =  5000000
    _shape = (2000, 3000)


class GmudGridTrishuli(_Grid):
    _p_str = """+proj=lcc +lat_1=30 +lat_2=62 +lat_0=0 +lon_0=105 +x_0=0 +y_0=0 
                +ellps=WGS84 +datum=WGS84 +units=m +no_defs"""
    _res= 1000
    _x0 = -2004000
    _y0 =  3877000

    _shape = (221, 168)

"""[(-2007500.0, 3881500.0),
    (-1847500.0, 3881500.0),
    (-1847500.0, 3677500.0),
    (-2007500.0, 3677500.0),
    (-2007500.0, 3881500.0)]
"""

class LatLonGrid(_Grid):
    _p_str = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    _res= 1
    _x0 = -1800
    _y0 = 900
    _shape = (1800, 3600)


_grids = {'gmuds':GmudsGrid, 'latlon':LatLonGrid, 'trishuli': GmudGridTrishuli}
get_grid = lambda x: _grids[x]


class _Dataset(object):

    data_vars = ['Tad', 'Pad', 'RHd', 'Ld']

    def __init__(self, nc_file, window=None, grid='gmuds', subset_grid='trishuli', load_on_init=False):
        self.grid = get_grid(grid)()
        self.nc_file = nc_file
        self.window = window
        nc_filename = os.path.basename(nc_file)
        nc_datestr = nc_filename.replace('.nc4', '')[-10:]
        self.start_date_time = datetime.datetime.strptime(nc_datestr, '%Y%m%d%H')
        self.subset_grid =  get_grid(subset_grid)()
        if load_on_init:
            self.data_array = self.load()

    def load(self):
        return self.subset_bygrid(self.subset_grid)

    def get_max_bounds(self, grid):
        x_, y_ = [], []
        for (x, y) in grid.coords:
            x_c, y_c = grid.transform(self.grid, x, y)
            x_.append(x_c)
            y_.append(y_c)
        return min(x_), max(x_), min(y_), max(y_)

    def subset_bygrid(self, grid):
        x0, x1, y0, y1 = self.get_max_bounds(grid)
        return self.subset(x0, x1, y0, y1)

    def read_dataset(self):

        ds = xr.open_dataset(self.nc_file)
        gt = (self.grid._x0, self.grid._res, 0, self.grid._y1, 0, -self.grid._res)
        try:
            self.missing_value = ds.MissingValue
            ds_ = {}
            # time treated differently in files, should remedy
            for var in self.data_vars:
                ds_[var] = (('y', 'x', 'time'), np.expand_dims(ds[var].T, -1))


            coords = {'easting' : (('x'),  ds['lon']),
                      'northing' : (('y'), ds['lat']),
                      'time': (('time'), np.array([self.start_date_time]).astype('datetime64[s]'))}
            data_set = xr.Dataset(ds_, coords=coords)
            data_set.attrs['gt'] = gt
            data_set.attrs['crs'] = self.grid._p_str
            for var in data_set.data_vars:
                data_set[var].attrs['gt'] = gt
                data_set[var].attrs['crs'] = self.grid._p_str
            return data_set
        except:
            raise()
        finally:
            ds.close()

    def subset(self, x0, x1, y0, y1, coords_grid=None):

        if coords_grid is not None:
            x0, y1 = coords_grid.transform(self.grid, x0, y1)
        col0, row0 = self.grid.inverse(x0, y1)

        if coords_grid is not None:
            x1, y0 = coords_grid.transform(self.grid, x1, y0)
        col1, row1 = self.grid.inverse(x1, y0)

        # gt for x, y on grid row, col coordinates,
        gt_x0, gt_y1 = self.grid.forward(col0, row0)
        gt = (gt_x0, self.grid._res, 0, gt_y1, 0, -self.grid._res)
        ds = xr.open_dataset(self.nc_file)

        #try:
        self.missing_value = ds.MissingValue
        ds_ = {}
        # time treated differently in files, should remedy
        for var in self.data_vars:
            ds_[var] = (('y', 'x', 'time'), np.expand_dims(ds[var][0,:,:].T[row0:row1, col0:col1], -1))
        coords = {'easting' : (('x'),  ds['lon'][col0:col1]),
                  'northing' : (('y'), ds['lat'][row0:row1]),
                  'time': (('time'), np.array([self.start_date_time]).astype('datetime64[s]'))}
        data_set = xr.Dataset(ds_, coords=coords)
        data_set.attrs['gt'] = gt
        data_set.attrs['crs'] = self.grid._p_str

        for var in data_set.data_vars:
            data_set[var].attrs['gt'] = gt
            data_set[var].attrs['crs'] = self.grid._p_str
        return data_set
        #except:
        #    raise()
        #finally:
        #    ds.close()

    def write_geotiff(self, data_array, file_name):
        assert data_array.ndim == 2
        _args = dict(driver='GTiff', height=data_array.shape[0], width=data_array.shape[1],
                     count=1, dtype=data_array.dtype.name, crs=data_array.attrs['crs'],
                     transform=Affine.from_gdal(*data_array.attrs['gt']))

        with rio.open(file_name, 'w', **_args) as geotif:
            geotif.write(data_array.values, 1)


_summary = {
    'all': ['Tad', 'Pad', 'RHd', 'Ld', 'Sd', 'U', 'V', 'PRECTOT'],
    'mean': ['Tad', 'Pad', 'RHd', 'Ld', 'Sd', 'U', 'V'],
    'sum': ['PRECTOT']
}

class Ambient(_Dataset):

    data_vars = ['Tad', 'Pad', 'RHd', 'Ld']


class ShortWave(_Dataset):

    data_vars = ['Sd']


class Wind(_Dataset):

    data_vars = ['U', 'V']


class Precipitation(_Dataset):

    data_vars = ['PRECTOT']


class DataLookup(object):

    # these lists define order
    data_names = ['ambient', 'shortwave', 'wind', 'precipitation']
    data_types = [Ambient, ShortWave, Wind, Precipitation]
    data_paths = [AMB_PATH, SWR_PATH, WND_PATH, PRC_PATH]

    def __init__(self, grid='gmuds'):
        for data_name, data_type, data_path in zip(self.data_names, self.data_types, self.data_paths):
            file_list = self.get_file_list(data_path)
            data_table = {}
            for nc_file in file_list:
                data_set = data_type(nc_file, load_on_init=False, grid=grid)
                data_table[data_set.start_date_time] = data_set
            setattr(self, data_name, data_table)

        data_collections, index = self.sort()

        self.data_collections = data_collections
        self.index = pd.DatetimeIndex(index)

        self.all_collections = copy.deepcopy(data_collections)
        self.all_index = copy.deepcopy(pd.DatetimeIndex(index))

    def sort(self):
        index = []
        for data_name in self.data_names:
            data_dates = list(getattr(self, data_name).keys())
            index.extend(data_dates)
        index = sorted(set(index))
        data_collections = []
        for date in index:
            data_collections.append(
                {
                    data_name: getattr(self, data_name)[date] for data_name in self.data_names
                }
            )
        return data_collections, index

    def get_file_list(self, data_path):
        file_list = []
        for root, folders, files in os.walk(data_path):
            for file in files:
                if file.endswith('.nc4'):
                    file_list.append(os.path.join(root, file))
        return file_list

    def groupby_date(self):
        seq = list(zip(self.data_collections, self.index.date))
        seq.sort(key=itemgetter(1))
        groups = groupby(seq, itemgetter(1))
        return [(key, [item[0] for item in data]) for (key, data) in groups]

    def select_hour(self, start_hour, stop_hour):
        date_index = np.where((self.index.hour >= start_hour) & (self.index.hour < stop_hour))[0]
        self.data_collections = self.data_collections[date_index]
        self.index = self.index[date_index]

    def select_date(self, date):
        date_index = np.where(self.index.date == date)[0]
        self.data_collections = self.data_collections[date_index]
        self.index = self.index[date_index]

    def reset(self):
        self.data_collections = copy.deepcopy(self.all_collections)
        self.index = copy.deepcopy(self.all_index)


class Merge(object):
    def __init__(self, data_collection):
        self.data_collection = data_collection
        # check if dates are equal
        self.dt = None
        for data in data_collection:
            if self.dt is None:
                self.dt = data_collection[data].start_date_time
            assert self.dt == data_collection[data].start_date_time
        self.data_set = self.combine()

    def combine(self):
        data_list = []
        for data in self.data_collection.values():
            print('Loading data: {}'.format(data.nc_file))
            data_list.append(data.load())
        #data_list = [data.load() for data in self.data_collection.values()]
        return xr.merge(data_list)


class TimeSeries(object):
    def __init__(self, merge_list):
        merge_data = [merge.data_set for merge in merge_list]
        self.data_set = xr.concat(merge_data, 'time')

    def write(self, file_name):
        self.data_set.to_netcdf(file_name)

    def write_daily(self, file_name):

        summary = {}
        summary['mean'] = self.data_set.mean(dim='time')
        summary['sum'] = self.data_set.sum(dim='time')

        daily_  = {}
        for var in _summary['all']:
            if var in _summary['mean']:
                daily_[var] = summary['mean'][var]

            elif var in _summary['sum']:
                daily_[var] = summary['sum'][var]

        daily_ds = xr.Dataset(daily_)
        daily_ds['PRECTOT'] = daily_ds.PRECTOT# * 86400 # to kg/m^2
        daily_ds.to_netcdf(file_name)


def write_daily():
    lu = DataLookup()
    group_dates = lu.groupby_date()
    for (date, collections) in group_dates:
        print('DATE: {}'.format(date))
        nc4_file =       'himat_gmuds_trishuli_hourly_{}.nc4'.format(date.strftime('%Y%m%d'))
        nc4_file_daily = 'himat_gmuds_trishuli_daily_{}.nc4'.format(date.strftime('%Y%m%d'))
        nc4_path = os.path.join(TRISHULI_PATH, nc4_file)
        nc4_path_daily = os.path.join(TRISHULI_PATH, nc4_file_daily)
        print(nc4_path)

        if not os.path.exists(nc4_path):
            collections_ = copy.deepcopy(collections)
            merge_list = [Merge(collection) for collection in collections_]
            time_series = TimeSeries(merge_list)
            time_series.write(nc4_path)
            time_series.write_daily(nc4_path_daily)


def write_latlon():

    # Trishuli Basin Lat/Lon Boundaries from shapefile (trishuli.shp)
    xMin, yMin = 84.9219, 27.7445
    xMax, yMax = 85.7988, 29.0059

    # Rounded lat/lon bounds (1/4 degree)
    # These will be different than actual bounds since they must
    # fall on GMUDS-GRID ROW.COL Grid (Azimuthal)
    xMin, yMin = 84.5, 27.5
    xMax, yMax = 86.0, 29.5
    tri_bounds_ll = [(xMin, yMin), (xMax, yMin), (xMax, yMax), (xMin, yMax), (xMin, yMin)]

    # get values from file
    nc_file = xr.open_dataset('./test/ts.nc4')
    easting_ = nc_file.easting.values
    northing_ = nc_file.northing.values
    gmuds_bounds = [(easting_[0], northing_[0]), (easting_[-1], northing_[0]),
                    (easting_[-1], northing_[-1]), (easting_[0], northing_[-1]),
                    (easting_[0], northing_[0])]

    # gmuds/latlon grid objects
    src_grid = get_grid('gmuds')()
    ll_grid =  get_grid('latlon')()
    dst_grid = get_grid('trishuli')()

    # calc lat/lon for gmud grid vertices
    easting, northing = np.meshgrid(easting_, northing_)
    lons, lats = ll_grid.itransform(src_grid, easting, northing)

    # GDAL geotransform

    gt_fromDat = (easting_[0], 1000, 0, northing_[0], 0, -1000)
    gt_fromObj = (dst_grid._x0, dst_grid._res, 0, dst_grid._y0, 0, -dst_grid._res)
    print('# to grid edges')
    print('geotransform = ({}, {}, {}, {}, {}, {})'.format(*gt_fromObj))
    print('# to grid centers')
    print('geotransform = ({}, {}, {}, {}, {}, {})'.format(*gt_fromDat))

    # to copy and paste into QGIS
    wkt_str_FromCalc = 'POLYGON(({}))'.format(','.join(['{} {}'.format(*ll_grid.itransform(src_grid, x, y)) for (x, y) in dst_grid.coords]))
    wkt_str_FromData = 'POLYGON(({}))'.format(','.join(['{} {}'.format(*ll_grid.itransform(src_grid, x, y)) for (x, y) in gmuds_bounds]))
    print(wkt_str_FromCalc)
    print(wkt_str_FromData)
    # write NetCDF
    latlon_ds = {}
    latlon_ds['latitude' ] = (('y', 'x'), lats)
    latlon_ds['longitude'] = (('y', 'x'), lons)
    coords = {'easting':  (('x'), easting_  ),
              'northing': (('y'), northing_ )}
    data_set = xr.Dataset(latlon_ds, coords=coords)
    data_set.to_netcdf('./test/latlon_trishuli_gmu_ds.nc4', mode='w')

    # write GeoTiff from obj
    _args = dict(driver='GTiff', height=dst_grid._shape[0], width=dst_grid._shape[1], count=2, dtype='float64',
                 crs=dst_grid._p_str,  transform=dst_grid.fwd)
    with rio.open('./test/test_geotiff_fromObj.tif', 'w', **_args) as geotif:
        geotif.write(lons, 1)
        geotif.write(lats, 2)

    # write GeoTiff from data
    _args = dict(driver='GTiff', height=northing_.size, width=easting_.size, count=2, dtype='float64',
                 crs=dst_grid._p_str,  transform=Affine.from_gdal(*gt_fromDat))
    with rio.open('./test/test_geotiff_fromDat.tif', 'w', **_args) as geotif:
        geotif.write(lons, 1)
        geotif.write(lats, 2)



"""
POLYGON((
        -1979779.8095711682 3710497.8493330223,
        -1895083.3404771092 3688871.543999353,
        -1860988.4414337012 3825577.9646937074,
        -1944161.1160320693 3846815.185979576,
        -1979779.8095711682 3710497.8493330223))
"""

def write_latlon_fromGrid(dst_name='gmuds', src_name='trishuli'):

    dst_grid = get_grid(dst_name)()
    src_grid = get_grid(src_name)()

    x0s, x1s, y0s, y1s = src_grid.get_max_bounds(dst_grid)

    col0, row0 = dst_grid.inverse(x0s, y1s)
    col1, row1 = dst_grid.inverse(x1s, y0s)

    # gt for x, y on grid row, col coordinates,
    gt_x0, gt_y1 = dst_grid.forward(col0, row0)

    gt = (gt_x0, dst_grid._res, 0, gt_y1, 0, -dst_grid._res)

    cols = np.arange(col0, col1, dtype='double')
    rows = np.arange(row0, row1, dtype='double')
    col_, row_ = np.meshgrid(cols, rows)
    col_ += .5
    row_ += .5

    lons, lats = dst_grid.forward(col_, row_, grid=LatLonGrid())

    easting, _ =  dst_grid.forward(cols, np.zeros_like(cols))
    _, northing = dst_grid.forward(np.zeros_like(rows), rows)

    latlon_ds = {}
    latlon_ds['latitude' ] = (('y', 'x'), lats)
    latlon_ds['longitude'] = (('y', 'x'), lons)
    coords = {'easting':  (('x'), easting ),
              'northing': (('y'), northing)}
    data_set = xr.Dataset(latlon_ds, coords=coords)
    data_set.to_netcdf(os.path.join(DAT_PATH, 'latlon_trishuli_gmu_ds.nc4'), mode='w')

def write_gtiff():

    nc_file = 'dat/gmu_trishuli/himat_gmuds_trishuli_dly_20071110.nc4'
    grid = GmudGridTrishuli()

    data_set = xr.open_dataset(nc_file)
    data_array = data_set.data_vars['Tad']

    _args = dict(driver='GTiff', height=data_array.shape[0], width=data_array.shape[1],
                 count=1, dtype=data_array.dtype.name, crs=grid.p.srs,
                 transform=Affine.from_gdal(*grid.gt))

    with rio.open(nc_file.replace('nc4', 'tif'), 'w', **_args) as geotif:
        geotif.write(data_array.values, 1)





WRITE = False # edit to change function

def main():
    if WRITE:
        write_daily()
    write_gtiff()

if __name__ == '__main__':
    sys.exit(main())
