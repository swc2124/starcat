from __future__ import division, absolute_import, print_function

import os
from time import sleep

import matplotlib.pyplot as plt
import ebf
import numpy as np

from c_functions import bin as _bin, integerize as _int

import printlib as plib
'''
'halo ebf directory': os.path.join(os.environ['HOMEPATH'], 'Desktop', 'halo_ebf'),
'halo filehandels': os.listdir(os.path.join(os.environ['HOMEPATH'], 'Desktop', 'halo_ebf')),
'''
program_data = {

    'halo file': 'not selected',
    'halo': 'not selected',
    'distance Mpc': 4.0,
    'app mag limit': 29.0,
    'cur dir': os.path.abspath(os.path.curdir),
    'processor count': os.environ['NUMBER_OF_PROCESSORS'],
    'regions': []

}


def find_data_dirs(pgrm_data=program_data):
    first = True
    for pth in os.path.abspath(os.path.curdir).split('\\')[:-1]:
        if first:
            data_dir = os.path.join(pth, os.path.sep)
            first = False
        data_dir = os.path.join(data_dir, pth)
        print(' 		', data_dir)

    data_dir = os.path.join(data_dir, 'data')
    print(' 		', data_dir)
    if not os.path.isdir(data_dir):
        os.mkdir(data_dir)
    pgrm_data['data dir'] = data_dir

    pgrm_data['grid dir'] = os.path.join(data_dir, 'grids')
    if not os.path.isdir(pgrm_data['grid dir']):
        os.mkdir(pgrm_data['grid dir'])

    pgrm_data['plot dir'] = os.path.join(data_dir, 'plots')
    if not os.path.isdir(pgrm_data['plot dir']):
        os.mkdir(pgrm_data['plot dir'])


# box coord set making function
def make_boxe(xy, size=100):
    '''

    '''
    x, y = xy
    x += 300
    y += 300
    segment = (size * 0.5)
    return ((x - segment, x + segment), (y - segment, y + segment))

# box placing function


def place_box(_ax, _box, _boxid):
    '''
    _ax is the subplot to plot on
    _box is ((x0, x1), (y0, y1))
    _boxid is the box ID

    make a single box from a set of numbers like this: box = ((x0, x1), (y0, y1))
    _ax.vlines(x0, y0, y1)
    _ax.vlines(x1, y0, y1)
    _ax.hlines(y0, x0, x1)
    _ax.hlines(y1, x0, x1)
    '''
    # unpack box tuples
    _age, _feh = _box

    # find box center
    a_center = _age[0] + ((_age[1] - _age[0]) * 0.5)
    f_center = _feh[1] + ((_feh[1] - _feh[0]) * 0.5)

    # set box label
    box_label = 'R ' + str(_boxid)
    _ax.text(a_center, f_center, box_label,
             ha='center',
             fontsize=10,
             color='k',
             bbox={'facecolor': 'white', 'alpha': 0.75, 'pad': 2})

    # params
    lwidth = 3
    lstyle = 'dashed'
    lclr = 'r'
    aph = .3

    # set lines
    _ax.vlines(_age[0], _feh[0], _feh[1],
               colors=lclr,
               linestyles=lstyle,
               alpha=aph,
               linewidth=lwidth)
    _ax.vlines(_age[1], _feh[0], _feh[1],
               colors=lclr,
               linestyles=lstyle,
               alpha=aph,
               linewidth=lwidth)
    _ax.hlines(_feh[0], _age[0], _age[1],
               colors=lclr,
               linestyles=lstyle,
               alpha=aph,
               linewidth=lwidth)
    _ax.hlines(_feh[1], _age[0], _age[1],
               colors=lclr,
               linestyles=lstyle,
               alpha=aph,
               linewidth=lwidth)


def fix_rslice(grid, rslices=[4]):

    # (1.0 * Config.getint('grid_options', 'size'))
    ratio = (2.0 * 300.0 / 600.0)
    x_center = (grid.shape[1] / 2)  # + (10 - x_bump[0])
    y_center = (grid.shape[0] / 2)  # + (10 - y_bump[0])
    for r in rslices:
        for i in range(grid.shape[0]):
            for q in range(grid.shape[1]):
                value = np.sqrt(
                    (np.square(i - y_center) + np.square(q - x_center)))
                value /= ratio
                if value > 300.0:
                    value = 0.0
                grid[i, q, r] = value
    return grid


def plot_halo(halo_fh, pgrm_data=program_data):
    '''
    halo = ebf.read(os.path.join(
        program_data['halo ebf directory'], program_data['halo file']))
    px, py = _int(halo['px'].astype(np.float64), halo['py'].astype(np.float64))
    ab_mags = halo['wfirst-hst_h158'].astype(np.float64)
    ap_mags = ab_mags + (5 * np.log10(4.0 * 1e5))
    r_proj = np.sqrt(np.square(halo['px']) + \
                     np.square(halo['py'])).astype(np.float64)
    lims = np.array([27.0, 28.0, 29.0], dtype=np.float64)
    grid = _bin(px, py, ab_mags, ap_mags, r_proj,
                lims, halo['satid'].astype(np.int32))
    '''
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.set_title(halo_fh)
    ticks = [str(i) for i in np.linspace(-150, 150, 7)]
    ax.set_xlabel('KPC')
    ax.set_xticklabels(ticks)
    ax.set_ylabel('KPC')
    ax.set_yticklabels(ticks)

    grid_fh = os.path.join(pgrm_data['grid dir'], halo_fh)
    print(' 		loading :', grid_fh)
    grid = np.load(grid_fh)
    print(' 		done')
    # fill radius gridslice
    print(' 		correcting radius')
    grid = fix_rslice(grid)
    print(' 		done')
    print(' 		making heat map')
    ax.pcolormesh(np.log10(grid[8:, :-3, 0]),
                  cmap=plt.cm.bone_r, vmin=1.0, vmax=4.5)
    print(' 		making contour plot')
    cp = ax.contour(grid[:, :, 4], [50, 100, 150, 300], colors='k',
                    linewidths=1.5, alpha=.25, linestyles='dashed')
    print(' 		setting contour labels')
    cl = ax.clabel(cp, [50, 100, 150, 300], inline=1, fmt='%s Kpc',
                   fontsize=10, color='k', linewidth=50, alpha=1)
    print(' 		setting limits')
    # limits
    center = grid.shape[0] / 2
    include = center / 2
    lim0 = center + include
    lim1 = center - include
    ax.set_xlim([lim1, lim0])
    ax.set_ylim([lim1, lim0])

    if pgrm_data['regions']:
        print(' 		setting boxes')
        # boxes [(((x1, x2), (y1, y2)), size)]
        for box in pgrm_data['regions']:
            print(' 		', box)
            place_box(ax, box, 1)

    print(' 		setting kpc grid')
    ax.axes.grid(alpha=.4, linestyle='dashed', color='grey')
    plot_fh = os.path.join(
        pgrm_data['plot dir'], halo_fh.split('_')[0] + '.png')
    print(' 		saving plot now to :', plot_fh)
    fig.savefig(plot_fh, dpi=800)
    program_data['plot'] = plot_fh
    print(' 		done')


def exit():
    selection = raw_input(' 		are you sure you wantto quit? y/[n]')
    if selection in ['y', 'Y', 'yes', 'Yes', 'YES']:
        import sys
        clear()
        sys.exit(0)
    main_menu()


def clear():
    # clear terminal
    if os.name == 'nt':
        clr_cmd = 'cls'
    else:
        clr_cmd = 'clear'
    os.system(clr_cmd)


def main_menu(pgrm_data=program_data):
    clear()
    print(plib.main_title)
    print(plib.main_menu_title)

    for key in pgrm_data.keys():
        if key == 'halo filehandels':
            continue
        space_to_add = 15 - len(key)
        print(' 		', key, ' ' * space_to_add, ' : ', pgrm_data[key])

    print('\n 		select from the following items')
    print(' 		----------------------------------------')
    for i, _menu in enumerate(all_menus):
        print(' 		[', i, '] ', _menu[0])
    selection = int(raw_input('\n 		menu number: '))

    clear()
    print(plib.main_title)
    all_menus[selection][1]()


def select_halo_file(pgrm_data=program_data):
    halo_grid_dir = program_data['grid dir']
    halo = program_data['halo']

    print('\n 		selecting halo file')
    print(' 		----------------------------------------')
    print(' 		current halo           :', halo)
    print(' 		current grid directory :', halo_grid_dir)

    halo_filehandels = os.listdir(halo_grid_dir)
    print('\n 		select halo file from current grid files')
    print(' 		----------------------------------------')
    for i, fh in enumerate(halo_filehandels):
        print(' 		[', i, '] ', fh.split('_')[0])
    selection = int(raw_input(' 		select halo: '))
    program_data['halo file'] = halo_filehandels[selection]
    halo_fh = program_data['halo file']
    program_data['halo'] = halo_fh.split('_')[0]
    halo = program_data['halo']
    print('\n 		----------------------------------------')
    print(' 		current halo           :', halo)
    print(' 		new halo file          :', halo_fh)
    print(' 		----------------------------------------')
    # TODO load saved regions
    pgrm_data['regions'] = []

    print('\n 		plotting halo now')
    plot_halo(halo_fh)
    main_menu()


def select_region(pgrm_data=program_data):
    clear()
    print(plib.main_title)
    print('\n 		selecting regions for ', pgrm_data['halo'])
    print(' 		----------------------------------------')
    print('\n 		current regions')
    print(' 		----------------------------------------')
    if not len(pgrm_data['regions']):
        print(' 		 none')
    else:
        for i, reg in enumerate(pgrm_data['regions']):
        	print(' 		[', i, '] ', reg)
    print('\n 		select from the following items')
    print(' 		----------------------------------------')
    for i, _menu in enumerate(select_region_menue):
        print(' 		[', i, '] ', _menu[0])
    selection = int(raw_input('\n 		menu number: '))
    select_region_menue[selection][1]()


def remove_region(pgrm_data=program_data):
    selection = int(raw_input('\n 		enter region number to delete'))
    del pgrm_data['regions'][selection]
    plot_halo(program_data['halo file'])
    select_region()


def new_region(pgrm_data=program_data):
    print('\n 		enter new region location and size')
    print(' 		----------------------------------------')
    new_x = int(raw_input(' 		new region X :'))
    print(' 		', new_x)
    new_y = int(raw_input(' 		new region Y :'))
    print(' 		', new_y)
    new_size = raw_input(' 		new region size : (enter for default)')
    if new_size == '':
        new_size = 10
    else:
        new_size = int(new_size)

    print(' 		', new_size)
    new_box = make_boxe((new_x, new_y), size=new_size)
    print(' 		new box:', new_box, 'size:', new_size)
    pgrm_data['regions'].append(new_box)
    plot_halo(program_data['halo file'])
    select_region()


all_menus = [
    ('select halo file', select_halo_file),
    ('select regions', select_region),
    ('exit', exit)
]

select_region_menue = [
	('remove_region', remove_region),
	('new_region', new_region),
	('main menu', main_menu)
]

if __name__ == '__main__':
    clear()
    for key in os.environ.keys():
        if key == 'PATH':
            print(key)
            for pth in os.environ[key].split(';'):
                print(' --> :', pth)
        else:
            space_to_add = 35 - len(key)
            print(key, ' ' * space_to_add, ' : ', os.environ[key])
    raw_input('\npress enter to continue')
    find_data_dirs()
    main_menu()
