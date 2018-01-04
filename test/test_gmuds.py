import os
import sys

import gmuds as gm

def test_read_nc():
    nc_file = 'test/Atmd2008051603.nc4'
    array = gm.Array(nc_file)
    print(array.nc_file)

def test_subset_array():
    nc_file = 'test/Atmd2008051603.nc4'
    gm.read_nc(nc_file)


if __name__ == '__main__':
    test_read_nc()