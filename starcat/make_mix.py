#! /usr/bin/env python
'''
Read in stellar catalogs from simulations with Absolute AB mags and
background galaxy catlogs. Using the WingsTips lib, produce mixed list
of objects including stars and appropriate sampling of background
galaxies in STIPS input format
'''

from wingtips import WingTips as stips
from wingtips import ascii
from wingtips import np
from wingtips import time

files = ['h15.shell.1Mpc.in', 'h15.shell.3Mpc.in',
             'h15.shell.5Mpc.in', 'h15.shell.10Mpc.in']
ZP_AB = np.array([26.365, 26.357, 26.320, 26.367, 25.913])
filters = ['Z087', 'Y106', 'J129', 'H158', 'F184']


def make_stips():
    '''
    Returns:
        NoneType
    '''
#   Step through files in "files".
    for i, infile in enumerate(files):

#       File handle prefix for later (starpre).
        starpre = '_'.join(infile.split('.')[:-1])

#       Extract distance Mpc from filename (dist).
        dist = float(infile.split('.')[2][:-3])

#       Read in the data file (data) and print.
        data = ascii.read(infile)
        print('\nRead in %s \n' % infile)

#       Right ascension (RA) & declination (DEC).
        RA = data['col1']
        DEC = data['col2']

#       Magnitude array (M).
        M = np.array([
            data['col3'], data['col4'],
            data['col5'], data['col6'],
            data['col7']]).T

#       What is this doing?
        temp = [
            stips.from_scratch(
                flux=stips.get_counts(
                    M[:, j],
                    ZP_AB[j],
                    dist=dist),
                ra=RA,
                dec=DEC,
                outfile=starpre + '_' + filt[0] + '.tbl')
            for j, filt in enumerate(filters)]
    return None


def mix_stips(_fltrs=filters, _fnames=files, _outprefix='Mixed'):
    '''[summary]



    Keyword Arguments:
        _fltrs {[type]}    (default: {filters})
        _fnames {[type]}   (default: {files})
        _outprefix {str}   (default: {'Mixed'})

    Returns:
        NoneType
    '''
#   Make an empty list to stash galaxies into (galaxies).
    galaxies = []

#   Step through files in "_fnames".
    for i, infile in enumerate(_fnames):

#       File handle prefix for later (starpre).
        starpre = '_'.join(infile.split('.')[:-1])

#       New empty list (radec).
        radec = []

#       Step through filters "_fltrs".
        for j, filt in enumerate(_fltrs):

            stars = stips([starpre + '_' + filt[0] + '.tbl'])

#           Is this the first infile?
            if i == 0:
                galaxies.append(stips([filt + '.txt']))
                galaxies[j].flux_to_Sb()

            if len(radec) == 0:
                radec = galaxies[j].random_radec_for(stars)

#
            galaxies[j].replace_radec(radec)

#
            stars.merge_with(galaxies[j])

#           # File handle (outfile).
            outfile = '_'.join((_outprefix, starpre, filt[0])) + '.tbl'

#           # Write stars to file.
            stars.write_stips(outfile, ipac=True)

#           Write new catalog header to file.
            with open(outfile, 'r+') as f:
                content = f.read()
                f.seek(0, 0)
                f.write(
                    '\\type = internal' + '\n' +
                    '\\filter = ' + str(filt) + '\n' +
                    '\\center = (' + str(stars.center[0]) +
                    '  ' + str(stars.center[1]) + ')\n' + content)
    return None

if __name__ == '__main__':
    tic = time.time()
    assert 3 / 2 == 1.5, 'Not running Python3 may lead to wrong results'
    make_stips()
    mix_stips()
    print('\n\nCompleted in %.3f seconds \n' % (time.time() - tic))
