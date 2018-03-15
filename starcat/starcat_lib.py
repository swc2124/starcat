# ============================================================================
# Author                : swc21
# Date                  : 2018-03-13 20:11:26
# Project               : GitHub
# File Name             : starcat_lib
# Last Modified by      : swc21
# Last Modified time    : 2018-03-14 12:26:30
# ============================================================================
# 

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import os

from astropy.table import Table

def select_stars(_stars, _region):
    # Find the stars in the table that
    # are in the region (idx).print('finding stars within region boundaries',
    # file=os.sys.stdout)
    print('x0:', _region['x0'], file=os.sys.stdout)
    print('x1:', _region['x1'], file=os.sys.stdout)
    print('y0:', _region['y0'], file=os.sys.stdout)
    print('y1:', _region['y1'], file=os.sys.stdout)
    idx_ = np.nonzero(
        np.logical_and(
            np.logical_and(
                _stars['x_int'] >= _region['x0'],
                _stars['x_int'] <= _region['x1']),
            np.logical_and(
                _stars['y_int'] >= _region['y0'],
                _stars['y_int'] <= _region['y1'])))[0]
    return _stars[idx_], idx_


def load_table(_table_fh):
    return Table.read(_table_fh, path='data')


def check_dir(_pth):
    if not os.path.isdir(_pth):
        os.mkdir(_pth)
        print(_pth, file=os.sys.stdout)


def my_mkdir(_dir):
    if not type(_dir) == list:
        _dir = [_dir]
    for _path in _dir:
        if not os.path.isdir(_path):
            print('making ', _path, file=os.sys.stderr)
            os.mkdir(_path)
    print('done', file=os.sys.stdout)


def catalog_path(_dir, _name, _r, _typ):
    print('making catalog directory', file=os.sys.stdout)
    out_dir = os.path.join(_dir, _name)
    check_dir(out_dir)
    out_dir = os.path.join(out_dir, _typ[0])
    check_dir(out_dir)
    return os.path.join(out_dir, _r + _typ[1])


class CatalogManager(Table):
    """docstring for CatalogManager"""

    def __init__(self, _output_dir=None):
        super(CatalogManager, self).__init__()
        self.catalogs = {}
        self.output_dir = _output_dir
        self.n_catalogs = 0
        self.n_regions = 0

    def new_catalog(self, _name):
        # _name would be like "halo02"
        self.catalogs[_name] = {'n_regions': 0}
        self.n_catalogs = len(self.catalogs)

    def new_region(self, _name, _region):
        # _region is a dict
        # _name is like "halo02"
        self.catalogs[_name][_region['name']] = _region
        self.catalogs[_name]['n_regions'] += 1

    def remove_region(self, _name, _region_name):
        del self.catalogs[_name][_region_name]
        self.catalogs[_name]['n_regions'] -= 1

    def make_all(self):
        return None

    @staticmethod
    def remove_all():
        return None




