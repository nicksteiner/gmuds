import os
import sys
import pickle

import numpy as np
import pyproj
import xarray as xr
from affine import Affine
import rasterio as rio

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
                self._y0, self._y0 +self.y_size)

    @property
    def coords(self):
        return [
            (self._x0, self._y0),
            (self._x0 + self.x_size, self._y0),
            (self._x0 + self.x_size, self._y0 + self.y_size),
            (self._x0, self._y0 + self.y_size),
            (self._x0, self._y0)
        ]


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
        return int(np.floor(col)), int(np.ceil(row))  #


class GmudsGrid(_Grid):
    _p_str = """+proj=lcc +lat_1=30 +lat_2=62 +lat_0=0 +lon_0=105 +x_0=0 +y_0=0 
                +ellps=WGS84 +datum=WGS84 +units=m +no_defs"""
    _res = 1000
    _x0 = -3000000
    _y0 = 5000000
    _shape = (2000, 3000)

class UtmGrid(_Grid):
    _p_str ="+proj=utm +zone=45 +datum=WGS84 +units=m +no_defs"
    _res= 30
    _x0 = 120000
    _y0 = 3300000
    _shape = (10000, 16000)

class UtmGridTrishuli(_Grid):
    _p_str ="+proj=utm +zone=45 +datum=WGS84 +units=m +no_defs"
    _res= 30
    _x0 = 270000
    _y0 = 3060000
    _shape = (6000, 4000)

class LatLonGrid(_Grid):
    _p_str = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    _res= 1
    _x0 = -1800
    _y0 = 900
    _shape = (1800, 3600)

_grids = {'gmuds':GmudsGrid, 'utm':UtmGrid, 'latlon':LatLonGrid}
get_grid = lambda x: _grids[x]

class Dataset(object):

    def __init__(self, nc_file, window=None, grid='gmuds'):
        self.grid = get_grid(grid)()
        self.nc_file = nc_file
        self.window = window

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

    def subset(self, x0, x1, y0, y1, coords_grid=None):

        if coords_grid is not None:
            x0, y1 = coords_grid.transform(self.grid, x0, y1)
        col0, row0 = self.grid.inverse(x0, y1)

        if coords_grid is not None:
            x1, y0 = coords_grid.transform(self.grid, x1, y0)
        col1, row1 = self.grid.inverse(x1, y0)

        gt = (x0, self.grid._res, 0, y1, 0, -self.grid._res)

        ds = xr.open_dataset(self.nc_file)

        try:
            ds_ = {}
            # time treated differently in files, should remedy
            for var in ds.data_vars:
                ds_[var] = (('y', 'x'), ds[var][0, row0:row1, col0:col1])  # row-col ordered
            coords = {'lon': (('x'), ds['lon'][col0:col1]), 'lat': (('y'), ds['lat'][row0:row1])}
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


    def write_geotiff(self, data_array, file_name):

        _args = dict(driver='GTiff', height=data_array.shape[0], width=data_array.shape[1],
                     count=1, dtype=data_array.dtype.name, crs=data_array.attrs['crs'],
                     transform=Affine.from_gdal(*data_array.attrs['gt']))

        with rio.open(file_name, 'w', **_args) as geotif:
            geotif.write(data_array.values, 1)


