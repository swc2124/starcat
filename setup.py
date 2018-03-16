
from numpy import get_include
from setuptools import setup

setup(
    name='starcat',
    version='0.1.0',
    author='Sol W. Courtney',
    author_email='swc2124@Columbia.edu',
    maintainer='Sol W. Courtney',
    maintainer_email='swc2124@Columbia.edu',
    url='https://github.com/swc2124/starcat',
    description=(
        'Python tools for working with stellar data '
        'in the form of numpy arrays.  Intended to '
        'be used for defining WFIRST survey strategies'
    ),
    download_url='https://github.com/swc2124/starcat.git',
    license='MIT',
    include_package_data=True,
    include_dirs=[get_include()],
    package_data={
        # If any package contains *.txt files, include them:
        '': ['*.hdf5', '.npy'],
        # And include any *.dat files found in the 'data' subdirectory
        # of the 'mypkg' package, also:
        'starcat': ['data/*'],
    },
    install_requires=[
        'docutils>=0.3',
        'numpy>=0.x',
        'ebfpy>=0.x',
        'astropy>=0.x',
        'matplotlib>=0.x'])
