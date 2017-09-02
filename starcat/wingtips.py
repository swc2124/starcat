#! /usr/bin/env python
'''WFIRST Infrared Nearby Galaxies Test Image Product Simulator.
Produces input files for the WFIRST STIPS simulator

[description]
'''
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os

from builtins import int
from io import open

import numpy as np
import time

from astropy import wcs
from astropy.io import ascii as astro_ascii
from astropy.io import fits
from astropy.table import Table

import starcat as _sc


class WingTips:
    '''
    WFIRST Infrared Nearby Galaxies Test Image Product Simulator
    Produces input files for the WFIRST STIPS simulator
    '''

    def __init__(self, infile=[], center=[0, 0]):
        '''[summary]

        [description]

        Parameters
        ----------
        0] : {[type]}
            [description]
        infile : {list}, optional
            [description] (the default is [], which [default_description])
        center : {list}, optional
            [description] (the default is [0, which [default_description])

        Returns
        -------
        [type]
            [description]
        '''
        if len(infile) == 0:
            self.tab = np.array([])
        else:
            if isinstance(infile, str):
                infile = [infile]
            self.tab = WingTips.read_stips(infile[0])
            if len(infile) > 1:
                for i in range(1, len(infile)):
                    _row = WingTips.read_stips(infile[i])
                    self.tab = np.vstack((self.tab, _row))
            center = WingTips.get_center(self.tab[:, 0], self.tab[:, 1])
        self.center = center
        self.n = self.tab.shape[0]
        self.infile = infile
        return None

    ''' Strip coordinates from WingTips object '''

    def strip_radec(self, hasID=False):
        '''[summary]

        [description]

        Parameters
        ----------
        hasID : {bool}, optional
            [description] (the default is False, which [default_description])

        Returns
        -------
        [type]
            [description]
        '''
        _i = int(hasID)
        self.tab = np.delete(self.tab, [_i, _i + 1], 1)
        return None

    ''' Attach given RA-DEC to WingTips object'''

    def attach_radec(self, radec, hasID=False):
        '''[summary]

        [description]

        Parameters
        ----------
        radec : {[type]}
            [description]
        hasID : {bool}, optional
            [description] (the default is False, which [default_description])

        Returns
        -------
        [type]
            [description]

        Raises
        ------
        ValueError
            [description]
        '''
        if self.n != radec.shape[0]:
            raise ValueError('Number of RA-DEC does not match sources')
        _i = int(hasID)
        self.tab = np.insert(self.tab, _i, radec.T, 1)
        self.center = WingTips.get_center(radec[:, 0 + _i], radec[:, 1 + _i])
        return None

    ''' Replace RA-DEC of WingTips object '''

    def replace_radec(self, radec, hasID=False):
        '''[summary]

        [description]

        Parameters
        ----------
        radec : {[type]}
            [description]
        hasID : {bool}, optional
            [description] (the default is False, which [default_description])

        Returns
        -------
        [type]
            [description]
        '''
        self.strip_radec(hasID)
        self.attach_radec(radec, hasID)
        return None

    '''
    Return random RA-DEC for given image or WingTips object
    Optionally, specify center and image size desired
    '''

    def random_radec_for(self, other, shape=(4096, 4096), sample=False, n=0, hasID=False):
        '''[summary]

        [description]

        Parameters
        ----------
        other : {[type]}
            [description]
        shape : {tuple}, optional
            [description] (the default is (4096, 4096), which [default_description])
        sample : {bool}, optional
            [description] (the default is False, which [default_description])
        n : {number}, optional
            [description] (the default is 0, which [default_description])
        hasID : {bool}, optional
            [description] (the default is False, which [default_description])

        Returns
        -------
        [type]
            [description]
        '''
        _i = int(hasID)
        try:
            if other.endswith('.fits'):
                return WingTips.random_radec(self.n, imfile=other)
        except AttributeError:
            if not sample:
                return WingTips.random_radec(self.n, center=other.center)
            elif not bool(n):
                return WingTips.sample_radec(n=self.n, radec1=False, radec2=other.tab[:, _i:_i + 1])
            else:
                return WingTips.sample_radec(n=n, radec1=self.tab[:, _i:_i + 1], radec2=other.tab[:, _i:_i + 1])

    ''' Merge two WingTips objects '''

    def merge_with(self, other, hasRADEC=True, hasID=False):
        '''[summary]

        [description]

        Parameters
        ----------
        other : {[type]}
            [description]
        hasRADEC : {bool}, optional
            [description] (the default is True, which [default_description])
        hasID : {bool}, optional
            [description] (the default is False, which [default_description])

        Returns
        -------
        [type]
            [description]

        Raises
        ------
        ValueError
            [description]
        '''
        if self.tab.shape[1] != other.tab.shape[1]:
            raise ValueError('Number of columns does not match',
                             self.tab.shape[1], other.tab.shape[1])
        self.tab = np.vstack((self.tab, other.tab))
        self.n = self.tab.shape[0]
        self.infile.append(other.infile)
        _i = int(hasID)
        if hasRADEC:
            self.center = WingTips.get_center(
                self.tab[:, 0 + _i], self.tab[:, 1 + _i])
        return None

    ''' Convert flux to surface brightness for sersic profile galaxies '''

    def flux_to_Sb(self, hasRADEC=True, hasID=False):
        '''[summary]

        [description]

        Parameters
        ----------
        hasRADEC : {bool}, optional
            [description] (the default is True, which [default_description])
        hasID : {bool}, optional
            [description] (the default is False, which [default_description])

        Returns
        -------
        [type]
            [description]
        '''
        _i = int(hasID)
        if hasRADEC:
            _i = _i + 2
        _f = self.tab[:, _i].astype(float)
        _r = self.tab[:, _i + 3].astype(float)
        _a = self.tab[:, _i + 5].astype(float)
        _s = (0.5 * _f) / (np.pi * _r**2 * _a)
        self.tab = np.delete(self.tab, _i, 1)
        self.tab = np.insert(self.tab, _i, _s.T, 1)
        return None

    ''' Write out a STIPS input file '''

    def write_stips(self, outfile='temp.txt', hasID=False, hasCmnt=False, saveID=False, ipac=False):
        '''[summary]

        [description]

        Parameters
        ----------
        outfile : {str}, optional
            [description] (the default is 'temp.txt', which [default_description])
        hasID : {bool}, optional
            [description] (the default is False, which [default_description])
        hasCmnt : {bool}, optional
            [description] (the default is False, which [default_description])
        saveID : {bool}, optional
            [description] (the default is False, which [default_description])
        ipac : {bool}, optional
            [description] (the default is False, which [default_description])

        Returns
        -------
        [type]
            [description]
        '''
        _tab = WingTips.get_tabular(self.tab, hasID, hasCmnt, saveID)
        _nms = ('id',  'ra',    'dec',   'flux',  'type',
                'n',   're',    'phi',   'ratio', 'notes')
        _fmt = ('%10d', '%15.7f', '%15.7f', '%15.7f', '%8s',
                '%10.3f', '%15.7f', '%15.7f', '%15.7f', '%8s')
        _t = Table(_tab, names=_nms)
        if ipac:
            astro_ascii.write(_t, outfile, format='ipac',
                              formats=dict(zip(_nms, _fmt)))
        else:
            astro_ascii.write(_t, outfile, format='fixed_width',
                              delimiter='', formats=dict(zip(_nms, _fmt)))
        return print('Wrote out %s \n' % outfile)

    ''' Build a WingTips class object from scratch '''
    @staticmethod
    def from_scratch(flux, ra=[], dec=[], center=[], ID=[], Type=[], n=[], re=[], phi=[], ratio=[], notes=[], outfile=''):
        '''[summary]

        [description]

        Parameters
        ----------
        flux : {[type]}
            [description]
        ra : {list}, optional
            [description] (the default is [], which [default_description])
        dec : {list}, optional
            [description] (the default is [], which [default_description])
        center : {list}, optional
            [description] (the default is [], which [default_description])
        ID : {list}, optional
            [description] (the default is [], which [default_description])
        Type : {list}, optional
            [description] (the default is [], which [default_description])
        n : {list}, optional
            [description] (the default is [], which [default_description])
        re : {list}, optional
            [description] (the default is [], which [default_description])
        phi : {list}, optional
            [description] (the default is [], which [default_description])
        ratio : {list}, optional
            [description] (the default is [], which [default_description])
        notes : {list}, optional
            [description] (the default is [], which [default_description])
        outfile : {str}, optional
            [description] (the default is '', which [default_description])

        Returns
        -------
        [type]
            [description]

        Raises
        ------
        ValueError
            [description]
        '''
        _temp = WingTips()
        _temp.n = len(flux)
        _temp.infile = ['fromScratch']

        if len(center) > 0:
            _temp.center = center
            if len(ra) == 0:
                radec = _temp.random_radec_for(_temp)
                ra, dec = radec[:, 0], radec[:, 1]
        elif ((len(ra) == len(dec)) & (len(ra) > 0)):
            _temp.center = WingTips.get_center(np.array(ra), np.array(dec))
        else:
            raise ValueError('Provide valid coordinate or center')

        if ((len(Type) == 0) | (Type is 'point') | (Type is 'sersic')):
            if ((len(Type) == 0) | (Type is 'point')):
                Type = np.repeat(np.array(['point']), len(flux))
                _ones = np.ones_like(flux)
                n, re, phi, ratio = _ones, _ones, _ones, _ones
            elif (Type == 'sersic'):
                Type = np.repeat(np.array(['sersic']), len(flux))
        elif (len(Type) == len(flux)):
            Type = np.array(Type)

        _tab = np.array([ra, dec, flux, Type, n, re, phi, ratio]).T

        if (len(ID) == len(flux)):
            _tab = np.hstack((np.array(ID, ndmin=2).T, _tab))
        if (len(notes) == len(flux)):
            _tab = np.hstack((_tab, np.array(notes, ndmin=2).T))

        _temp.tab = np.array(_tab)

        if outfile is '':
            return _temp
        else:
            _temp.write_stips(outfile, hasID=bool(
                ID), hasCmnt=bool(notes), saveID=bool(ID))
            return None

    '''
    Read in a STIPS input file in astro_ascii format and
    return corrsponding NumPy array
    '''
    @staticmethod
    def read_stips(infile, getRADEC=True, getID=False, getCmnt=False):
        '''[summary]

        [description]

        Parameters
        ----------
        infile : {[type]}
            [description]
        getRADEC : {bool}, optional
            [description] (the default is True, which [default_description])
        getID : {bool}, optional
            [description] (the default is False, which [default_description])
        getCmnt : {bool}, optional
            [description] (the default is False, which [default_description])

        Returns
        -------
        [type]
            [description]
        '''
        _tab = []
        _infile = astro_ascii.read(infile)
        print('\nRead in %s \n' % infile)

        if getID:
            _tab.append(_infile['id'])
        if getRADEC:
            _tab.append(_infile['ra'])
            _tab.append(_infile['dec'])

        _tab.append(_infile['flux'])
        _tab.append(_infile['type'])
        _tab.append(_infile['n'])
        _tab.append(_infile['re'])
        _tab.append(_infile['phi'])
        _tab.append(_infile['ratio'])

        if getCmnt:
            _tab.append(_infile['comment'])

        return np.array(_tab).T

    ''' Return tabular lists for STIPS input file columns '''
    @staticmethod
    def get_tabular(_tab, hasID=False, hasCmnt=False, saveID=False):
        '''[summary]

        [description]

        Parameters
        ----------
        _tab : {[type]}
            [description]
        hasID : {bool}, optional
            [description] (the default is False, which [default_description])
        hasCmnt : {bool}, optional
            [description] (the default is False, which [default_description])
        saveID : {bool}, optional
            [description] (the default is False, which [default_description])
        '''
        _i = int(hasID)
        if ~saveID:
            _n = _tab.shape[0]
            _ID = np.array(np.linspace(1, _n, _n), ndmin=2).T
            _tab = np.hstack((_ID, _tab[:, _i:]))
        if ~hasCmnt:
            _cmnt = np.array(
                np.repeat(np.array(['comment']), _tab.shape[0],), ndmin=2).T
            _tab = np.hstack((_tab, _cmnt))
        return [_tab[:, 0].astype(float), _tab[:, 1].astype(float), _tab[:, 2].astype(float),
                _tab[:, 3].astype(float), _tab[:, 4], _tab[:, 5].astype(float),
                _tab[:, 6].astype(float), _tab[:, 7].astype(float),
                _tab[:, 8].astype(float), _tab[:, 9]]

    ''' Build WCS coordinate system from scratch '''
    @staticmethod
    def create_wcs(centers=[0, 0], crpix=[2048, 2048], cdelt=[-0.11 / 3600, 0.11 / 3600], cunit=['deg', 'deg'],
                   ctype=['RA---TAN', 'DEC--TAN'], lonpole=180, latpole=24.333335,
                   equinox=2000.0, radesys='ICRS'):
        '''[summary]

        [description]

        Parameters
        ----------
        0] : {[type]}
            [description]
        2048] : {[type]}
            [description]
        0.11 / 3600] : {[type]}
            [description]
        'deg'] : {[type]}
            [description]
        'DEC--TAN'] : {[type]}
            [description]
        centers : {list}, optional
            [description] (the default is [0, which [default_description])
        crpix : {list}, optional
            [description] (the default is [2048, which [default_description])
        cdelt : {list}, optional
            [description] (the default is [-0.11 / 3600, which [default_description])
        cunit : {list}, optional
            [description] (the default is ['deg', which [default_description])
        ctype : {list}, optional
            [description] (the default is ['RA---TAN', which [default_description])
        lonpole : {number}, optional
            [description] (the default is 180, which [default_description])
        latpole : {number}, optional
            [description] (the default is 24.333335, which [default_description])
        equinox : {number}, optional
            [description] (the default is 2000.0, which [default_description])
        radesys : {str}, optional
            [description] (the default is 'ICRS', which [default_description])

        Returns
        -------
        [type]
            [description]
        '''
        _w = wcs.WCS()
        _w.wcs.cdelt = cdelt
        _w.wcs.crpix = crpix
        _w.wcs.crval = centers
        _w.wcs.cunit = cunit
        _w.wcs.ctype = ctype
        _w.wcs.lonpole = lonpole
        _w.wcs.latpole = latpole
        _w.wcs.radesys = radesys
        _w.wcs.equinox = equinox
        return _w

    ''' Return coordinate system for given image file'''
    @staticmethod
    def read_wcs(imfile):
        '''[summary]

        [description]

        Parameters
        ----------
        imfile : {[type]}
            [description]

        Returns
        -------
        [type]
            [description]
        '''
        print('Getting coordinates from %s \n' % imfile)
        return wcs.WCS(fits.open(imfile)[1].header)

    ''' Return 'n' random radec for given image file or coordinate list '''
    @staticmethod
    def random_radec(n=10, center=[0, 0], shape=(4096, 4096), imfile=''):
        '''[summary]

        [description]

        Parameters
        ----------
        0] : {[type]}
            [description]
        n : {number}, optional
            [description] (the default is 10, which [default_description])
        center : {list}, optional
            [description] (the default is [0, which [default_description])
        shape : {tuple}, optional
            [description] (the default is (4096, 4096), which [default_description])
        imfile : {str}, optional
            [description] (the default is '', which [default_description])

        Returns
        -------
        [type]
            [description]
        '''
        _xy = np.random.rand(n, 2) * shape
        if imfile is not '':
            _w = WingTips.read_wcs(imfile)
        else:
            _w = WingTips.create_wcs(center)
        return _w.wcs_pix2world(_xy, 1)

    '''
    Return a random sample of 'n' RA-DEC coordinates from 'radec2'
    If radec1 is specified, then replace 'n' radom coordinates
    in 'radec1' with random sample from 'radec2'
    '''
    @staticmethod
    def sample_radec(n=10, radec1=False, radec2=[]):
        '''[summary]

        [description]

        Parameters
        ----------
        n : {number}, optional
            [description] (the default is 10, which [default_description])
        radec1 : {bool}, optional
            [description] (the default is False, which [default_description])
        radec2 : {list}, optional
            [description] (the default is [], which [default_description])

        Returns
        -------
        [type]
            [description]
        '''
        in2 = np.random.randint(0, radec2.shape[0], n)
        if ~radec1:
            return radec2[in2, :]
        else:
            in1 = np.random.randint(0, radec1.shape[0], n)
            radec1[in1, :] = radec2[in2, :]
            return radec1

    ''' Return mean of RA-DEC positions given '''
    @staticmethod
    def get_center(ra, dec):
        '''[summary]

        [description]

        Parameters
        ----------
        ra : {[type]}
            [description]
        dec : {[type]}
            [description]
        '''
        return [ra.astype(float).mean(), dec.astype(float).mean()]

    '''
    Convert mags to WFI instrument counts
    Default is apparent AB mags
    Specify 'dist' if absolute mags
    Specify AB_Vega if Vega mags
    '''
    @staticmethod
    def get_counts(mag, ZP, dist=0, AB_Vega=0):
        '''[summary]

        [description]

        Parameters
        ----------
        mag : {[type]}
            [description]
        ZP : {[type]}
            [description]
        dist : {number}, optional
            [description] (the default is 0, which [default_description])
        AB_Vega : {number}, optional
            [description] (the default is 0, which [default_description])

        Returns
        -------
        number
            [description]
        '''
        if bool(dist):
            print('\nDistance is d = %4.2f Mpc\n' % dist)
            u = 25 + 5 * np.log10(dist)
            mag = mag + u
        if bool(AB_Vega):
            mag = mag + AB_Vega
        return 10**((mag - ZP) / (-2.5))
