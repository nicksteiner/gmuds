import os
import sys
import pickle

import numpy as np
import pyproj
import xarray as xr
from affine import Affine

# class
class _Grid(object):
    def __init__(self):
        self.p = pyproj.Proj(self._p_str)
        self.gt = (self._x0, self._res, 0, self._y0, 0, -self._res)
        self.fwd = Affine.from_gdal(*self.gt)
        self.rev = ~self.fwd

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
        return int(col), int(row)


class GmudsGrid(_Grid):
    _p_str = """+proj=lcc +lat_1=30 +lat_2=62 +lat_0=0 +lon_0=105 +x_0=0 +y_0=0 
                +ellps=WGS84 +datum=WGS84 +units=m +no_defs"""
    _res = 1000
    _x0 = -3000000
    _y0 = 5000000
    _size = (2000, 3000)

class UtmGrid(_Grid):
    _p_str ="+proj=utm +zone=45 +datum=WGS84 +units=m +no_defs"
    _res= 30
    _x0 = 120000
    _y0 = 3300000
    _size = (10000, 16000)

class LatLonGrid(_Grid):
    _p_str = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    _res= 1
    _x0 = -1800
    _y0 = 900
    _size = (1800, 3600)


class Array(object):

    grid = GmudsGrid()

    def __init__(self, nc_file, window=None):
        self.nc_file = nc_file
        self.window = window

    def subset(self, x0, x1, y0, y1, grid=None):
        ds = xr.open_dataset(self.nc_file)
        try:
            ds
        finally:
            ds.close()