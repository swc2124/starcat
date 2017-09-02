#
# -*- coding: utf-8 -*-
# GitHub  setup
# =============================================================================
# Name                : setup.py
# Date                : 2017-08-19 18:54:50
# Author              : sol courtney
# GitHub              : https://github.com/swc2124
# Affiliation         : Columbia University NYC, NY
# Email               : swc2124@columbia.edu
# Language            : Python
# Last Modified by    : swc21
# Last Modified time  : 2017-08-31 09:44:37
# =============================================================================

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os

from numpy import get_include
from setuptools import find_packages
from setuptools import setup


def list_data_files():
    _top_dir = os.path.join(os.path.curdir, 'starcat') + os.path.sep
    _include_list = []
    _data_dirs = ['tables', 'grids']
    _file_types = ['.txt', '.hdf5', '.in', '.ebf', '.npy']
    _pkg_data_path = os.path.join(os.path.curdir, 'starcat', 'data')
    for directory in _data_dirs:
        for name in os.listdir(os.path.join(_pkg_data_path, directory)):
            lcl_pth = os.path.join(_pkg_data_path, directory, name)
            if os.path.isdir(lcl_pth):
                for name_ in os.listdir(lcl_pth):
                    lcl_lcl_pth = os.path.join(lcl_pth, name_)
                    if os.path.isfile(lcl_lcl_pth):
                        if os.path.splitext(lcl_lcl_pth)[-1].lower() in _file_types:
                            _include_list.append(
                                lcl_lcl_pth.replace(_top_dir, ''))
            else:
                if os.path.splitext(lcl_pth)[-1].lower() in _file_types:
                    _include_list.append(lcl_pth.replace(_top_dir, ''))
    return _include_list

setup(

    name='starcat',
    version='0.1.0',

    author='Sol W. Courtney',
    author_email='swc2124@Columbia.edu',

    maintainer='Sol W. Courtney',
    maintainer_email='swc2124@Columbia.edu',

    url='https://github.com/swc2124/starcat',

    description=(
        'Python tools for working with stellar data.'
        'Intended to be used for defining WFIRST survey strategies'
    ),
    download_url='https://github.com/swc2124/starcat.git',
    license='MIT',
    include_package_data=True,
    include_dirs=[get_include()],
    packages=['starcat'],
    package_dir={'starcat': 'starcat'},
    package_data={'starcat': list_data_files()},
    python_requires='>=2.6, !=3.0.*, !=3.1.*, !=3.2.*, <4',
    classifiers=[

        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: WFIRST :: Simulation',

        # Pick your license as you wish (should match "license" above)
        'License :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    install_requires=[
        'pillow>=4.0.0',
        'future>=0.16.0',
        'six>=1.1.0',
        'docutils>=0.13.1',
        'numpy>=1.11.2',
        'ebfpy>=0.0.14',
        'astropy>=1.3.2',
        'matplotlib>=2.0.2'],
    entry_points={
        'console_scripts': [
            'starcat = starcat.starcat_gui:run_gui', ]
    }
)
