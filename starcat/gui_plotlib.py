
from __future__ import division
from __future__ import print_function

import matplotlib as mpl
import numpy as np
import os
mpl.use('TkAgg')

import matplotlib.pyplot as plt

from astropy.table import Table


def kpc_to_degree(d_mpc):
    d = 1e3  # Kpc
    D = d_mpc * 1e6
    degree_mod = 57.295  # 180.0 / np.pi
    return degree_mod * (d / D)


def kpc_to_arcmin(d_mpc):
    d = 1e3  # Kpc
    D = d_mpc * 1e6
    arcmin_mod = 3437.746  # (60.0 * 180.0) / np.pi
    return arcmin_mod * (d / D)


def kpc_to_arcsec(d_mpc):
    d = 1e3  # Kpc
    D = d_mpc * 1e6
    arcsec_mod = 206264.806  # (3600.0 * 180.0) / np.pi
    return arcsec_mod * (d / D)


def make_box(session, name=None):
    '''
    box maker

    function to help the region selection function.
    makes a  set of box coord for the plot function to use.

    Arguments:
        session.current_x  center x coord for new box
        session.current_y  center y coord for new box
        session.region_size  ize of box in grid units

    '''
    #
    if str(session.halo) in session.regions.keys():
        session.regions[str(session.halo)][0] += 1
        session.regions[str(session.halo)][1] += 1
    else:
        session.regions[str(session.halo)] = [0, 0]

    # ratio of arcsec/Kpc
    arcsec_per_kpc_ratio = kpc_to_arcsec(session.distance)
    print('arcsec_per_kpc_ratio:', arcsec_per_kpc_ratio)

    # the inverse is Kpc/arcsec
    kpc_per_arcsec_ratio = 1.0 / arcsec_per_kpc_ratio
    print('kpc_per_arcsec_ratio:', kpc_per_arcsec_ratio)

    # this is how many x,y pixels there are
    xpix, ypix = session.region_ccd_size
    print('xpix, ypix: ', xpix, ypix)

    # to get arcsecs we multiply Npixel (xpix) * 0.11 arcsec/pixel
    # (region_pixel_size)
    detector_arcsec_x = xpix * session.region_pixel_size
    detector_arcsec_y = ypix * session.region_pixel_size
    print('detector_arcsec_x:', detector_arcsec_x)
    print('detector_arcsec_y:', detector_arcsec_y)

    # to get Kpc for the grid we divide by kpc_to_arcsec()
    # so its arcsec * Kpc/arcsec
    detector_kpc_x = detector_arcsec_x * kpc_per_arcsec_ratio
    detector_kpc_y = detector_arcsec_y * kpc_per_arcsec_ratio
    print('detector_kpc_x:', detector_kpc_x)
    print('detector_kpc_y:', detector_kpc_y)

    # half distance
    segment_x = detector_kpc_x * 0.5
    segment_y = detector_kpc_y * 0.5
    print('segment_x:', segment_x)
    print('segment_y:', segment_y)

    # current x, y locations of the mouse click (center of the detector)
    x_now = session.current_x
    y_now = session.current_y
    print('x_now:', x_now)
    print('y_now:', y_now)

    # (x0, x1) & (y0, y1)
    # x and y are backward in tkinter
    # that's why the signs are opposite
    xbox = (x_now - segment_x, x_now + segment_x)
    ybox = (y_now + segment_y, y_now - segment_y)
    print('xbox:', xbox)
    print('ybox:', ybox)

    # make a new region dict and append it to this halo's list of regions
    new_region = {}
    new_region['halo'] = str(session.halo)
    if name:
        new_region['name'] = name + str(session.regions[str(session.halo)][0])
    else:
        new_region['name'] = 'R ' + str(session.regions[str(session.halo)][0])

    new_region['box'] = (xbox, ybox)

    # need to offset for the grid
    # the grid is 600 x 600 so we add 1/2 i.e. 300 (grid_center_x,
    # grid_center_y)
    new_region['x0'] = xbox[0] + session.grid_center_x
    new_region['x1'] = xbox[1] + session.grid_center_x
    new_region['y0'] = ybox[1] + session.grid_center_y
    new_region['y1'] = ybox[0] + session.grid_center_y

    # append this new region dict to this halo's list of regions
    session.regions[str(session.halo)].append(new_region)


def fix_rslice(grid, d_mpc=4.0, unit='kpc', rslices=[4]):
    '''
    helper function for plot_halo()

    fills the entire grid slice for radius
    this is needed for a smooth contour over plot

    Arguments:
        grid {np.ndarray} -- the grid

    Keyword Arguments:
        rslices {list} -- slice with radius for all stars (default: {[14]})

    Returns:
        np.ndarray -- the same grid but with filled radial slice
    '''
    # TODO - determine the ratio
    # ratio = (2.0 * 300.0 / 600.0)

    # find center of grid
    x_center = (grid.shape[1] / 2)
    y_center = (grid.shape[0] / 2)

    # conver to units
    if unit == 'arcmin':
        ratio = kpc_to_arcmin(d_mpc)
    elif unit == 'arcsec':
        ratio = kpc_to_arcsec(d_mpc)
    elif unit == 'degree':
        ratio = kpc_to_degree(d_mpc)
    else:
        ratio = 1

    # iterate over whole grid one by one
    for r in rslices:
        for i in range(grid.shape[0]):
            for q in range(grid.shape[1]):
                value = np.sqrt((
                    np.square(i - y_center) +
                    np.square(q - x_center))
                )
                if value > 300.0:
                    value = 0.0
                value *= ratio
                grid[i, q, r] = value

    return grid

    # end of function =====================================================


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
    print('----------------------')
    print('place_box()')
    # unpack box tuples
    print('box:', _box)
    xlims, ylims = _box
    print('xlims, ylims:', xlims, ylims)

    # find box center
    x_center = xlims[0] + ((xlims[1] - xlims[0]) * 0.5)
    y_center = ylims[1] - ((ylims[1] - ylims[0]) * 0.5)

    print('x_center:', x_center)
    print('y_center:', y_center)

    # set box label
    if type(_boxid) == int:
        box_label = 'R ' + str(_boxid)
        print('box label:', box_label)
        _ax.text(
            x_center,
            y_center,
            box_label,
            ha='center',
            fontsize=10,
            color='k',
            bbox={'facecolor': 'white',
                  'alpha': 0.75,
                  'pad': 2
                  }
        )
    print('----------------------')
    # box params
    lwidth = 5
    lstyle = 'dashed'
    lclr = 'r'
    aph = 1

    # set lines
    _ax.vlines(xlims[0], ylims[0], ylims[1],
               colors=lclr,
               linestyles=lstyle,
               alpha=aph,
               linewidth=lwidth)
    _ax.vlines(xlims[1], ylims[0], ylims[1],
               colors=lclr,
               linestyles=lstyle,
               alpha=aph,
               linewidth=lwidth)
    _ax.hlines(ylims[0], xlims[0], xlims[1],
               colors=lclr,
               linestyles=lstyle,
               alpha=aph,
               linewidth=lwidth)
    _ax.hlines(ylims[1], xlims[0], xlims[1],
               colors=lclr,
               linestyles=lstyle,
               alpha=aph,
               linewidth=lwidth)
    return _ax

    # end of function =====================================================


def save_colorbar(_vmin, _vmax, _label, _fh, _cmap=mpl.cm.bone_r):
    fig = plt.figure(figsize=(8, 1))
    ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
    norm = mpl.colors.Normalize(vmin=_vmin, vmax=_vmax)
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=_cmap,
                                    norm=norm,
                                    orientation='horizontal')
    cb1.set_label(_label)
    fig.savefig(_fh)


def plot_cmd(session):
    '''
    plot a cmd of the local region
    '''
    make_box(session, 'CMD')

    region = session.regions[session.halo][-1]

    # make & set figure size
    fig = plt.figure(figsize=(10, 10))

    # make a single subplot
    ax = fig.add_subplot(111)

    # set titles & labels for axes & plot
    ax_title = str(session.halo) + ' ' + str(session.distance) + ' Mpc ' + \
        session.filter + '\nCMD' + str(session.regions[str(session.halo)][0])
    ax.set_title(ax_title)
    ax.set_xlabel('Z087 - H158', fontsize=7)
    ax.set_ylabel('H158', fontsize=7)
    idx = np.nonzero(
        np.logical_and(
            np.logical_and(
                session.table['x_int'] >= region['x0'],
                session.table['x_int'] <= region['x1']),
            np.logical_and(
                session.table['y_int'] >= region['y0'],
                session.table['y_int'] <= region['y1'])
        )
    )[0]

    mag_table = Table.read(session.mag_table_fh, format='hdf5', path='data')

    clr_arr = mag_table['m_z087'][idx] - mag_table['m_h158'][idx]
    clr_y_arr = mag_table['m_h158'][idx]

    ax.scatter(clr_arr, clr_y_arr, s=1, alpha=.75)

    ax.set_xlim([clr_arr.min(), clr_arr.max()])
    ax.set_ylim([clr_y_arr.max(), clr_y_arr.min()])

    # plot grid
    print('setting ' + session.plot_units + ' grid')
    ax.axes.grid(alpha=.4, linestyle='dashed', color='grey')

    # return fig

    # make plot filehandel (fh) and save plot
    plot_fh = os.path.join(session.plot_dir, session.halo +
                           region['name'] + '_cmd' + session.plot_extention)
    print('saving plot now to :', plot_fh)
    fig.savefig(plot_fh, dpi=session.plot_dpi)
    print('done')


def plot_halo(session):
    '''
    plot_halo is the only plot function used by starcat
    plots a single png file to specified directory
    '''
    # make & set figure size
    fig = plt.figure(figsize=(5, 5))

    # make a single subplot
    ax = fig.add_subplot(111)

    # set titles & labels for axes & plot
    ax_title = str(session.halo) + ' ' + \
        str(session.distance) + ' Mpc ' + session.filter
    ax.set_title(ax_title)
    ax.set_xlabel(session.plot_units, fontsize=7)
    ax.set_ylabel(session.plot_units, fontsize=7)

    # make & set tick labels for xy axes
    if session.plot_units == 'kpc':
        _levels = None
        clabel_frmt = '%s Kpc'
        fsize = 10
        lwidth = 45
        c_lwidth = 1.25
        mod = 1.0
    elif session.plot_units == 'arcmin':
        _levels = None
        clabel_frmt = '%s arcmin'
        fsize = 10
        lwidth = 45
        c_lwidth = .8
        mod = np.square(kpc_to_arcmin(session.distance))
    elif session.plot_units == 'arcsec':
        _levels = None
        clabel_frmt = '%s arcsec'
        fsize = 10
        lwidth = 45
        c_lwidth = .8
        mod = np.square(kpc_to_arcsec(session.distance))
    elif session.plot_units == 'degree':
        _levels = None
        clabel_frmt = '%s deg'
        fsize = 10
        lwidth = 45
        c_lwidth = 1.25
        mod = np.square(kpc_to_degree(session.distance))
    # make grid filehandel (fh) & load grid
    # fill radius grid slice for contour
    print('correcting radius in ', session.plot_units, 'units')
    print('distance Mpc:', session.distance)
    grid = fix_rslice(session.grid, session.distance, session.plot_units)

    # plot heat map
    plt_grid = np.log10(grid[8:, :-3, 0] / mod)
    min_max = np.unique(plt_grid.flatten())
    min_max.sort()
    _vmin = min_max[2]
    _vmax = min_max[-3]
    _cmap = plt.cm.bone_r
    print('making heat map')
    hm = ax.pcolormesh(plt_grid,
                       cmap=_cmap,
                       vmin=_vmin,
                       vmax=_vmax)

    # plot contour map to overlay on heat map
    print('making contour plot')
    #_levels = [50, 100, 150, 300]
    cp = ax.contour(grid[:, :, 4],
                    levels=_levels,
                    colors='k',
                    linewidths=c_lwidth,
                    alpha=.25,
                    linestyles='dashed')

    # labels for contour overlay
    print('setting contour labels')
    cl = ax.clabel(cp,
                   levels=_levels,
                   inline=1,
                   fmt='%s ' + session.plot_units,
                   fontsize=fsize,
                   color='k',
                   linewidth=lwidth,
                   alpha=1)

    # set plot xy limits
    print('setting limits')
    center = grid.shape[0] / 2
    include = session.plot_radius_kpc
    lim0 = center + include
    lim1 = center - include
    ax.set_xlim([lim1, lim0])
    ax.set_ylim([lim1, lim0])

    # plot region boxes for each region
    if str(session.halo) in session.regions.keys():
        print('setting boxes')
        for region in session.regions[str(session.halo)][2:]:

            # print(box)
            #ax = place_box(ax, box, i)
            print('------------------')
            print('place_box()')

            # unpack box tuples
            print('box:', region['box'])
            xlims, ylims = region['box']
            print('xlims, ylims:', xlims, ylims)

            # find box center
            x_center = xlims[0] + ((xlims[1] - xlims[0]) * 0.5)
            y_center = ylims[1] - 4.0 * ((ylims[1] - ylims[0]) * 0.5)

            print('x_center:', x_center)
            print('y_center:', y_center)

            # set box label
            box_label = region['name']
            print('box label:', box_label)
            ax.text(
                x_center + 300,
                y_center + 300,
                box_label,
                ha='center',
                fontsize=7,
                color='k',
                bbox={'facecolor': 'white',
                      'alpha': 0.75,
                      'pad': 2
                      }
            )
            print('------------------')
            # box params
            lwidth = 3
            lstyle = 'dashed'
            lclr = 'r'
            aph = .3

            x0 = region['x0']
            x1 = region['x1']
            y0 = region['y0']
            y1 = region['y1']

            # set lines
            ax.vlines(x0, y0, y1,
                      colors=lclr,
                      linestyles=lstyle,
                      alpha=aph,
                      linewidth=lwidth)
            ax.vlines(x1, y0, y1,
                      colors=lclr,
                      linestyles=lstyle,
                      alpha=aph,
                      linewidth=lwidth)
            ax.hlines(y0, x0, x1,
                      colors=lclr,
                      linestyles=lstyle,
                      alpha=aph,
                      linewidth=lwidth)
            ax.hlines(y1, x0, x1,
                      colors=lclr,
                      linestyles=lstyle,
                      alpha=aph,
                      linewidth=lwidth)

    # set unit tick labels
    current_ticks = ax.axes.get_xticks()
    c_ticks = current_ticks - center
    if session.plot_units == 'kpc':
        ticks = c_ticks

    elif session.plot_units == 'arcmin':
        ticks = (c_ticks * kpc_to_arcmin(session.distance)).round(1)

    elif session.plot_units == 'arcsec':
        ticks = (c_ticks * kpc_to_arcsec(session.distance)).round(1)

    elif session.plot_units == 'degree':
        ticks = (c_ticks * kpc_to_degree(session.distance)).round(1)

    ax.set_xticklabels(ticks, fontsize=10)
    ax.set_yticklabels(ticks, fontsize=10)

    # plot grid
    print('setting ' + session.plot_units + ' grid')
    ax.axes.grid(alpha=.4, linestyle='dashed', color='grey')

    # return fig

    # make plot filehandel (fh) and save plot
    plot_fh = session.plot_fh
    print('saving plot now to :', plot_fh)
    fig.savefig(plot_fh, dpi=session.plot_dpi)
    _label = 'Log(nstars / ' + session.plot_units + '^2)'
    _fh = os.path.join(session.plot_dir, session.halo +
                       '_cb' + session.plot_extention)
    session.plot_cb_fh = _fh
    session.plot_has_cb = True
    save_colorbar(_vmin, _vmax, _label, _fh)
    session.has_plot = True
    print('done')

    # end of function ========================================================
