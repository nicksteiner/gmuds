import os
import sys
import pickle

main_path = os.path.dirname(os.path.realpath(__file__)) + '/..'
sys.path.append(main_path)

import gmuds as gm



def test_gmuds_grid():
    gmuds_grid = gm.GmudsGrid()
    assert True


def test_gmuds_transform():
    gmuds_grid = gm.GmudsGrid()
    x, y = 100, 100
    x_, y_ = gmuds_grid.itransform(*[gmuds_grid] + list(gmuds_grid.transform(gmuds_grid, x, y)))
    assert gm.np.isclose((x_, y_), (x, y)).all()


def test_gmuds_forward():
    gmuds_grid = gm.GmudsGrid()
    col, row = 1500, 1000
    col_, row_ = gmuds_grid.inverse(*gmuds_grid.forward(col, row))
    assert gm.np.isclose((col_, row_), (col, row)).all()


def test_latlon2gmuds_transform():
    gmuds_grid = gm.GmudsGrid()
    latln_grid = gm.LatLonGrid()
    x, y = -1500000, 4000000
    lon, lat = latln_grid.itransform(gmuds_grid, x, y)
    x_, y_ = latln_grid.transform(gmuds_grid, lon, lat)
    assert gm.np.isclose((x_, y_), (x, y)).all()


def test_latlon2gmuds_inverse():
    gmuds_grid = gm.GmudsGrid()
    latln_grid = gm.LatLonGrid()
    #using transform
    lat, lon = 28.4044168199, 85.3545670252
    x, y = latln_grid.transform(gmuds_grid, lon, lat)
    row, col =  gmuds_grid.inverse(x, y)
    # using inverse
    row_, col_ = gmuds_grid.inverse(lon, lat, grid=latln_grid)
    assert gm.np.isclose((col_, row_), (col, row)).all()


def test_latlon2gmuds_forward():
    gmuds_grid = gm.GmudsGrid()
    latln_grid = gm.LatLonGrid()
    col, row = 1500, 1000
    #using transform
    x, y = gmuds_grid.forward(col, row)
    lon, lat = latln_grid.itransform(gmuds_grid, x, y)
    # using forward
    lon_, lat_ = gmuds_grid.forward(col, row, grid=latln_grid)
    assert gm.np.isclose((lon_, lat_), (lon, lat)).all()


def test_get_bounds():
    utm_grid = gm.UtmGrid()
    gmuds_grid = gm.GmudsGrid()
    latln_grid = gm.LatLonGrid()
    utm_grid_trishuli = gm.UtmGridTrishuli()


    utm_bnds = 270000, 390000, 3060000, 3240000
    utm_bounds = [
        (utm_bnds[0], utm_bnds[2]),
        (utm_bnds[0], utm_bnds[3]),
        (utm_bnds[1], utm_bnds[3]),
        (utm_bnds[1], utm_bnds[2]),
        (utm_bnds[0], utm_bnds[2])
    ]
    #print(utm_grid_trishuli.coords)

    #print(utm_bounds)

    latlon_bounds = [latln_grid.itransform(utm_grid, x, y) for x, y in utm_bounds]
    wkt_poly = 'POLYGON(({}))'.format(','.join(['{} {}'.format(*ln_lt) for ln_lt in latlon_bounds]))
    #print(wkt_poly)

    colrow_bounds = [gmuds_grid.inverse(x, y, grid=utm_grid) for x, y in utm_bounds]
    colrow_str  = 'Bounds: {}'.format(','.join(['{} {}'.format(*ln_lt) for ln_lt in colrow_bounds]))
    #print(colrow_str)

    utm_bounds_ = [gmuds_grid.forward(col, row, grid=utm_grid) for col, row in colrow_bounds]
    assert gm.np.isclose(utm_bounds_, utm_bounds, rtol=gmuds_grid._res/2).all()


def test_array():
    nc_file = 'test/Atmd2008051603.nc4'
    try:
        assert os.path.exists(nc_file)
    except:
        raise("Test data not loaded, too big for github: {}".format(nc_file))
    gm.Ambient(nc_file)
    assert True


def test_dataset_subset():
    nc_file = 'test/Atmd2008051603.nc4'
    try:
        assert os.path.exists(nc_file)
    except:
        raise("Test data not loaded, too big for github: {}".format(nc_file))
    array = gm.Ambient(nc_file)

    utm_grid_trishuli = gm.UtmGridTrishuli()
    x0, x1, y0, y1 = utm_grid_trishuli.bounds
    subset = array.subset(x0, x1, y0, y1, coords_grid=utm_grid_trishuli)
    assert True


def test_write_geotiff():
    nc_file = 'test/Atmd2008051603.nc4'
    d_var = 'Tad'
    try:
        assert os.path.exists(nc_file)
    except:
        raise("Test data not loaded, too big for github: {}".format(nc_file))
    array = gm.Ambient(nc_file)

    utm_grid_trishuli = gm.GmudGridTrishuli()
    x0, x1, y0, y1 = array.get_max_bounds(utm_grid_trishuli)
    subset = array.subset(x0, x1, y0, y1)
    #subset_ = array.subset_bygrid(utm_grid_trishuli)
    #assert (subset_['Tad'].values == subset['Tad'].values).all()
    data_array = subset[d_var][:,:,0]
    array.write_geotiff(data_array, 'test/test_{}.tif'.format(d_var))

    assert True

def test_write_full_geotiff():
    nc_file = 'test/Atmd2008051603.nc4'
    try:
        assert os.path.exists(nc_file)
    except:
        raise ("Test data not loaded, too big for github: {}".format(nc_file))
    dat_ = gm.xr.open_dataset(nc_file)
    _args = dict(driver='GTiff', height=2000, width=3000, count=1, dtype='float64',
                 crs='+proj=lcc +lat_1=30 +lat_2=62 +lat_0=0 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs',
                 transform=gm.Affine.from_gdal(*(-3000000, 1000, 0, 5000000, 0, -1000)))
    with gm.rio.open('./test/test_geotiff_full2.tif', 'w', **_args) as geotif:
        geotif.write(dat_.Tad.values[0, :, :].T, 1)
    assert True


def test_sw():
    nc_file = './test/SW2008012101.nc4'
    sw = gm.ShortWave(nc_file, load_on_init=True)


def test_wind():
    nc_file = './test/Wind2008012101.nc4'
    wn = gm.Wind(nc_file, load_on_init=True)


def test_merge():
    nc_file = './test/SW2008012101.nc4'
    sw = gm.ShortWave(nc_file, load_on_init=True)
    nc_file = './test/Wind2008012101.nc4'
    wn = gm.Wind(nc_file, load_on_init=True)
    nc_file = './test/Atmd2008012101.nc4'
    am = gm.Ambient(nc_file, load_on_init=True)
    _data_collection = {
        'ambient': am,
        'shortwave': sw,
        'wind': wn
    }

    full = gm.Merge(_data_collection)


def test_data_lookup():
    try:
        lu = gm.DataLookup()
        print(lu.index)
    except:
        print('Please run on server.')
        pass


def test_time_series():
    _data_collection = {
        'ambient': gm.Ambient('./test/Atmd2008012101.nc4', load_on_init=True),
        'shortwave': gm.ShortWave( './test/SW2008012101.nc4', load_on_init=True),
        'wind': gm.Wind('./test/Wind2008012101.nc4', load_on_init=True)
    }
    merge1 = gm.Merge(_data_collection)
    _data_collection = {
        'ambient': gm.Ambient('./test/Atmd2008012102.nc4', load_on_init=True),
        'shortwave': gm.ShortWave( './test/SW2008012102.nc4', load_on_init=True),
        'wind': gm.Wind('./test/Wind2008012102.nc4', load_on_init=True)
    }
    merge2 = gm.Merge(_data_collection)
    merge_list = gm.TimeSeries([merge1, merge2])
    pass


def test_lu_sort():
    pck_file = os.path.join(main_path, './test/lu.pck')
    lu = pickle.load(open(pck_file, 'rb'))
    group_dates = lu.groupby_date()
    print(lu.index)


def test_lu_merge():
    pck_file = os.path.join(main_path, './test/lu.pck')
    lu = pickle.load(open(pck_file, 'rb'))
    group_dates = lu.groupby_date()
    date_ = group_dates[0][0]
    merges_ = [gm.Merge(collection) for collection in group_dates[0][1]]
    assert True


def test_lu_time_series():
    pck_file = os.path.join(main_path, './test/lu.pck')
    lu = pickle.load(open(pck_file, 'rb'))
    group_dates = lu.groupby_date()
    date_ = group_dates[0][0]
    merges_ = [gm.Merge(collection) for collection in group_dates[0][1]]
    time_series = gm.TimeSeries(merges_)
    print('Time Series')
    pck_file = os.path.join(main_path, './test/ts.pck')
    print(pck_file)
    pickle.dump(time_series, open(pck_file, 'wb'))
    assert True


def write_lu_pickle():
    lu = gm.DataLookup()
    assert(len(lu.index) > 0)
    pck_file = os.path.join(main_path, './test/lu.pck')
    print(pck_file)
    pickle.dump(lu, open(pck_file, 'wb'))


def test_write_time_series():
    print('Time Series')
    pck_file = os.path.join(main_path, './test/ts.pck')
    print(pck_file)
    time_series = pickle.load(open(pck_file, 'rb'))
    nc4_file = os.path.join(main_path, './test/ts.nc4')

    time_series.write(nc4_file)


def test_write_latlon():
    gm.write_latlon()

if __name__ == '__main__':

    test_gmuds_grid()
    test_gmuds_transform()
    test_gmuds_forward()
    test_latlon2gmuds_transform()
    test_latlon2gmuds_inverse()
    test_latlon2gmuds_forward()
    test_get_bounds()
    test_array()
    test_dataset_subset()
    test_write_geotiff()
    test_sw()
    test_wind()
    test_merge()
    #test_data_lookup()
    #write_lu_pickle()
    #test_time_series()
    #test_get_daily()
    #test_lu_merge()
    #test_lu_time_series()
    #test_write_time_series()

    #test_write_latlon()

    pass