import os
import sys

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

    utm_bnds = 270000, 390000, 3060000, 3240000
    utm_bounds = [
        (utm_bnds[0], utm_bnds[2]),
        (utm_bnds[0], utm_bnds[3]),
        (utm_bnds[1], utm_bnds[3]),
        (utm_bnds[1], utm_bnds[2]),
        (utm_bnds[0], utm_bnds[2])
    ]
    #print(utm_bounds_)

    latlon_bounds = [latln_grid.itransform(utm_grid, x, y) for x, y in utm_bounds]
    wkt_poly = 'POLYGON(({}))'.format(','.join(['{} {}'.format(*ln_lt) for ln_lt in latlon_bounds]))
    #print(wkt_poly)

    colrow_bounds = [gmuds_grid.inverse(x, y, grid=utm_grid) for x, y in utm_bounds]
    colrow_str  = 'Bounds: {}'.format(','.join(['{} {}'.format(*ln_lt) for ln_lt in colrow_bounds]))
    #print(colrow_str)

    utm_bounds_ = [gmuds_grid.forward(col, row, grid=utm_grid) for col, row in colrow_bounds]
    assert gm.np.isclose(utm_bounds_, utm_bounds, rtol=gmuds_grid._res/2)


def test_array():
    pass


def test_array_subset():
    pass

if __name__ == '__main__':
    test_gmuds_grid()
    test_gmuds_transform()
    test_gmuds_forward()
    test_latlon2gmuds_transform()
    test_latlon2gmuds_inverse()
    test_latlon2gmuds_forward()
    test_get_bounds()