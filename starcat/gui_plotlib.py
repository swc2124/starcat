
from __future__ import division, print_function

import os
import numpy as np
import matplotlib
matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

import matplotlib.pyplot as plt

def kpc_to_degree(d_mpc):
    d = 1.0
    D = d_mpc * 1e6
    degree_mod = 180.0 / np.pi
    return degree_mod * (d / D)

def kpc_to_arcmin(d_mpc):
    d = 1.0
    D = d_mpc * 1e6
    arcmin_mod = (60.0 * 180.0) / np.pi
    return arcmin_mod * (d / D)

def kpc_to_arcsec(d_mpc):
    d = 1.0
    D = d_mpc * 1e6
    arcsec_mod = 206264.80624709636 #(3600.0 * 180.0) / np.pi
    return arcsec_mod * (d / D)

def make_box(session):
    '''
    box maker

    function to help the region selection function.
    makes a  set of box coord for the plot function to use.

    Arguments:
        session.current_x  center x coord for new box
        session.current_y  center y coord for new box
        session.region_size  ize of box in grid units

    '''
    # split the length
    segment = session.region_size * 0.5
    xbox = (session.current_x - segment, session.current_x + segment)
    ybox = (session.current_y + segment, session.current_y - segment)
    session.number_of_regions += 1
    session.number_of_regions_tot += 1
    new_region = {}
    new_region['halo'] = session.halo
    new_region['name'] = 'R ' + str(session.number_of_regions_tot)
    new_region['box'] = (xbox, ybox)
    new_region['x0'] = xbox[0] + 300
    new_region['x1'] = xbox[1] + 300
    new_region['y0'] = ybox[1] + 300
    new_region['y1'] = ybox[0] + 300
    session.regions.append(new_region)
    


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

    print('x_center:',x_center)
    print('y_center:',y_center)
    

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
    ax_title = session.halo + ' ' + str(session.distance) + ' Mpc ' + session.filter
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
    elif session.plot_units == 'arcmin':
        _levels = None
        clabel_frmt = '%s arcmin'
        fsize = 10
        lwidth = 45
        c_lwidth = .8
    elif session.plot_units == 'arcsec':
        _levels = None
        clabel_frmt = '%s arcsec'
        fsize = 10
        lwidth = 45
        c_lwidth = .8
    elif session.plot_units == 'degree':
        _levels = None
        clabel_frmt = '%s deg'
        fsize = 10
        lwidth = 45
        c_lwidth = 1.25
    # make grid filehandel (fh) & load grid
    # fill radius grid slice for contour
    print('correcting radius in ', session.plot_units, 'units')
    grid = fix_rslice(session.grid, session.distance, session.plot_units)

    # plot heat map
    print('making heat map')
    ax.pcolormesh(np.log10(grid[8:, :-3, 0]),
                  cmap=plt.cm.bone_r,
                  vmin=1.0,
                  vmax=4.5)

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
    if session.regions:
        print('setting boxes')
        for region in session.regions:
            #print(box)
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
        ticks = (c_ticks * kpc_to_arcmin(session.distance)).round(2)

    elif session.plot_units == 'arcsec':
        ticks = (c_ticks * kpc_to_arcsec(session.distance)).round(2)
 
    elif session.plot_units == 'degree':
        ticks = (c_ticks * kpc_to_degree(session.distance)).round(4)

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
    session.has_plot = True
    print('done')

    # end of function ========================================================
