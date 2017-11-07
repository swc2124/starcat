#
# -*- coding: utf-8 -*-
# GitHub  starcat_gui
# =============================================================================
# Name                : starcat_gui.py
# Date                : 2017-08-16 22:35:03
# Author              : sol courtney
# GitHub              : https://github.com/swc2124
# Affiliation         : Columbia University NYC, NY
# Email               : swc2124@columbia.edu
# Language            : Python
# Last Modified by    : swc21
# Last Modified time  : 2017-08-31 17:35:02
# =============================================================================
#
'''This is the main file to be run as python starcat_gui.py

Attributes
----------
root : {[type]}
    [description]
app : {[type]}
    [description]
root.mainloop() : {[type]}
    [description]

'''
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


import shutil

from builtins import dict
from builtins import int
from future.utils import iteritems
from io import open
from starcat_lib import os

try:
    from PIL import Image as PImage
    from PIL import ImageTk
    PhotoImage = ImageTk.PhotoImage
except ImportError as e:
    print(e[0], file=os.sys.stderr)
    from tkinter import PhotoImage
    from tkinter import Image as PImage

print(os.sys.version, file=os.sys.stdout)
if (os.sys.version_info > (3, 0)):
    # Python 3 code in this block
    import tkinter.filedialog as tkFileDialog
    import tkinter.simpledialog as tkSimpleDialog
else:
    # Python 2 code in this block
    import tkFileDialog
    import tkSimpleDialog

from tkinter import Frame
from tkinter import Label
from tkinter import Menu
from tkinter import PanedWindow
from tkinter import SUNKEN
from tkinter import Tk
from tkinter import W



from astropy.table import hstack
from starcat_lib import Table
from starcat_lib import np

try:
    from . import gui_plotlib as gplt
except ImportError as e:
    print(e[0], file=os.sys.stderr)
    import gui_plotlib as gplt


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


def output_catalog(_region, _c_typ, _out_pth):
    '''[summary]

    [description]

    Parameters
    ----------
    _region : {dict}
        dictionary of all the regions data.

    _c_typ : {list}
        _c_typ[0] = 'name'
        _c_typ[1] = '.ext'
        _c_typ[2] = 'format'

    _out_pth : {str}
        full PATH to output directory.
    '''
    # Make the file handle for the new catalog (cat_fh).
    cat_fh_ = catalog_path(_out_pth, _region['halo'], _region['name'], _c_typ)
    # Load halo data table (table)
    table_ = load_table(_region['table_fh'])
    # Find stars in region (region_idx).
    # Select stars from table (stars).
    stars_, r_idx_ = select_stars(table_, _region)
    # Make catalog.
    if _c_typ[0] == 'STIPS':
        pass
    else:
        # Load position table (pos_table).
        pos_table_ = load_table(_region['pxpy_table_fh'])
        # Load magnitude table (mag_table).
        mag_table_ = load_table(_region['mag_table_fh'])
        # Hstack the two together.
        fits_table_ = hstack([pos_table_[r_idx_], mag_table_[r_idx_], stars_])
        fits_table_.pprint()
        if 'spinbin_output_fh' in fits_table_.meta.keys():
            del fits_table_.meta['spinbin_output_fh']
        print('writing catalog to disc')
        fits_table_.write(cat_fh_, format=_c_typ[2], overwrite=True)
        print('adding row to record table')
        row_ = [_region['halo']]
        row_.append(_region['name'])
        row_.append(len(idx))
        row_.append(fits_table_['mact'].sum())
        row_.append(_region['x0'])
        row_.append(_region['x1'])
        row_.append(_region['y0'])
        row_.append(_region['y1'])
        row_.append(fits_table_['feh'].min())
        row_.append(fits_table_['feh'].mean())
        row_.append(fits_table_['feh'].max())
        row_.append(fits_table_['age'].min())
        row_.append(fits_table_['age'].mean())
        row_.append(fits_table_['age'].max())
        print(_region['name'], 'done', file=os.sys.stdout)
    return row_


class StatusBar(Frame):

    def __init__(self, master):
        Frame.__init__(self, master)
        self.label = Label(self, bd=1, relief=SUNKEN, anchor=W)
        self.label.grid()

    def set(self, format, *args):
        self.label.config(text=format % args)
        self.label.update_idletasks()

    def clear(self):
        self.label.config(text="")
        self.label.update_idletasks()


class App:

    def __init__(self, master):

        self.master = master
        self.master.title('StarCat')
        self.master.resizable(True, True)

        # file tree
        try:
            import starcat as _sc
            pkg_pth = os.path.dirname(_sc.__file__).lower()
            # preexisting directories
            self.data_dir = os.path.join(pkg_pth, 'data')
            self.table_dir = os.path.join(self.data_dir, 'tables')
            self.grid_dir = os.path.join(self.data_dir, 'grids')
            # New folders for output
            py_pth = os.sys.prefix.lower()
            cur_dir = os.path.abspath(os.path.curdir).lower()
            commonprefix = os.path.commonprefix([py_pth, pkg_pth, cur_dir])
            self.prefix = os.path.join(commonprefix, 'Starcat')
            self.save_dir = os.path.join(self.prefix, 'saved')
            self.plot_dir = os.path.join(self.prefix, 'plots')
            self.catalog_dir = os.path.join(self.prefix, 'catalogs')
            new_dirs = [self.prefix, self.plot_dir, self.catalog_dir]
            my_mkdir(new_dirs)

        except (ImportError, IOError, WindowsError) as e:
            #print(e[0], file=os.sys.stderr)
            os.sys.exit(e)

        #
        self.filter_type = 'h158'
        # regions
        self.regions = {}
        # 0.11 arcsec WFI detector pixel scale in arcseconds/pixel
        self.region_pixel_size = 0.11
        # (x, y) Number of active detector columns
        self.region_ccd_size = (4088, 4088)

        # status messages
        self.smsg_region = '[current regions: 0]    '
        self.smsg_mouse_xy = '[] '
        self.smsg_kpc_xy = '[]   '

        # plot display
        self.plot_size = (800, 800)
        self.plotwindow_size = 1000
        self.has_plot = False
        self.show_cb = False
        min_w = self.plot_size[0] + 50
        min_h = self.plot_size[0] + 100
        master.minsize(width=min_w, height=min_h)

        # plot settings
        self.plot_dpi = 400
        self.plot_units = 'kpc'
        self.plot_known_units = ['kpc', 'arcmin', 'arcsec', 'degree']
        self.plot_extention = '.png'
        self.plot_radius_kpc = 150
        self.plot_has_cb = False

        # catalog settings
        self.catalog_known_ext = ['.fits', '.hdf5', '.csv', '.txt']
        self.catalog_extention = '.fits'

        # status bar
        self.statusbar = StatusBar(self.master)
        self.statusbar.grid(row=2, column=0, columnspan=2)

        # plot window
        self.make_plotwindow()

        # top menu
        topmenu = Menu(self.master)
        topmenu.config(bd=5, title='Main Menu')
        self.master.config(menu=topmenu)

        # top menu - file menu
        file_menu = Menu(topmenu)
        file_menu.add_command(
            label="Exit",
            command=self.master.quit)
        topmenu.add_cascade(
            label="File",
            menu=file_menu)

        # top menu - halo menu
        halo_menu = Menu(topmenu)
        halo_menu.add_command(
            label="Select Halo",
            command=self.load_halo)
        topmenu.add_cascade(
            label="Halo",
            menu=halo_menu)

        # top menu - plot menu
        plot_menu = Menu(topmenu)
        plot_menu.add_command(
            label="Replot",
            command=self.plot_halo)
        plot_menu.add_separator()
        plot_menu.add_command(
            label="Set plot directory",
            command=self.plot_set_directory)

        # top menu - plot menu - dpi
        plot_dpi_menu = Menu(plot_menu)
        plot_dpi_menu.add_command(
            label="100",
            command=self.plot_set_dpi_100)
        plot_dpi_menu.add_command(
            label="200",
            command=self.plot_set_dpi_200)
        plot_dpi_menu.add_command(
            label="400",
            command=self.plot_set_dpi_400)
        plot_dpi_menu.add_command(
            label="600",
            command=self.plot_set_dpi_600)
        plot_dpi_menu.add_command(
            label="800",
            command=self.plot_set_dpi_800)
        plot_menu.add_cascade(
            label='Set dpi',
            menu=plot_dpi_menu)

        # top menu - plot menu - units
        plot_units_menu = Menu(plot_menu)
        plot_units_menu.add_command(
            label='Kpc',
            command=self.plot_set_units_kpc)
        plot_units_menu.add_command(
            label='Arcsec',
            command=self.plot_set_units_arcsec)
        plot_units_menu.add_command(
            label='Arcmin',
            command=self.plot_set_units_arcmin)
        plot_units_menu.add_command(
            label='Degree',
            command=self.plot_set_units_deg)
        plot_menu.add_cascade(
            label='Set Units',
            menu=plot_units_menu)

        # top menu - plot menu - radius
        plot_radius_menu = Menu(plot_menu)
        plot_radius_menu.add_command(
            label="150",
            command=self.plot_set_radius_150)
        plot_radius_menu.add_command(
            label="200",
            command=self.plot_set_radius_200)
        plot_radius_menu.add_command(
            label="250",
            command=self.plot_set_radius_250)
        plot_radius_menu.add_command(
            label="300",
            command=self.plot_set_radius_300)
        plot_menu.add_cascade(
            label='Set Radius Kpc',
            menu=plot_radius_menu)

        # top menu - plot menu - extension
        plot_ext_menu = Menu(plot_menu)
        plot_ext_menu.add_command(
            label="png",
            command=self.plot_set_ext_png)
        plot_ext_menu.add_command(
            label="pdf",
            command=self.plot_set_ext_pdf)
        plot_ext_menu.add_command(
            label="jpeg",
            command=self.plot_set_ext_jpeg)
        plot_menu.add_cascade(
            label='Set Extension',
            menu=plot_ext_menu)

        plot_menu.add_separator()
        plot_menu.add_command(
            label='Clear plot',
            command=self.clear_label_image,
            background='orange2')
        plot_menu.add_separator()
        plot_menu.add_command(
            label="Remove all",
            command=self.plot_remove_all,
            background='OrangeRed2')
        topmenu.add_cascade(
            label="Plot",
            menu=plot_menu)

        # top menu - region menu
        region_menu = Menu(topmenu)
        detector_menu = Menu(region_menu)
        detector_menu.add_command(
            label='H4RG-10',
            command=self.h4rg_10)
        region_menu.add_cascade(
            label='Detectors',
            menu=detector_menu)
        region_menu.add_separator()
        region_menu.add_command(
            label='Detector Size',
            command=self.region_set_ccd_size)
        region_menu.add_command(
            label='Pixel Size',
            command=self.region_set_pix_size)
        topmenu.add_cascade(
            label='Region',
            menu=region_menu)

        # top menu - catalog menu
        self.catalog_menu = Menu(topmenu)
        self.catalog_menu.add_command(
            label='Make Catalogs',
            command=self.catalog_makeall)
        catalog_type_menu = Menu(catalog_menu)
        catalog_type_menu.add_command(
            label='STIPS',
            command=self.catalog_set_STIPS)
        catalog_type_menu.add_command(
            label='fits',
            command=self.catalog_set_fits)
        catalog_menu.add_cascade(
            label='Format',
            menu=catalog_type_menu)

        self.catalog_menu.add_command(
            label='Set catalog output directory',
            command=self.catalog_set_output_dir)

        self.catalog_menu.add_separator()
        self.catalog_menu.add_command(
            label='Remove all',
            command=self.catalog_remove_all,
            background='OrangeRed2')
        topmenu.add_cascade(
            label='Catalog',
            menu=self.catalog_menu)

        # right click menu
        self.rightclick = Menu(self.plotwindow, tearoff=1)
        self.rightclick.add_command(
            label='Select Region',
            command=self.place_region)
        self.rightclick.add_command(
            label='Remove region',
            command=self.region_remove)
        self.rightclick.add_separator()
        self.rightclick.add_command(
            label='Make CMD',
            command=self.make_cmd)
        self.rightclick.add_separator()
        self.rightclick.add_command(
            label='Show Color-bar',
            command=self.plot_set_show_colorbar,
            background='lightgreen')

    def record_table(self):
        _names = ['halo', 'region', 'n_stars', 'M_tot',
                  'x0', 'x1', 'y0', 'y1',
                  'feh_min', 'feh_mean', 'feh_max',
                  'age_min', 'age_mean', 'age_max']
        _dtype = ['S10', 'S10', 'i', 'f',
                  'i', 'i', 'i', 'i',
                  'f', 'f', 'f',
                  'f', 'f', 'f']
        return Table(names=_names, dtype=_dtype)

    def catalog_set_output_dir(self):
        self.catalog_dir = tkFileDialog.askdirectory(
            initialdir=os.path.curdir,
            mustexist=True,
            parent=self.plotwindow,
            title='Catalog Output Directory')

    def catalog_remove_all(self):
        for file in os.listdir(self.catalog_dir):
            print(' --> evaluating', file)
            if file[:4] == 'halo':
                dir_path = os.path.join(self.catalog_dir, file)
                if os.path.isdir(dir_path):
                    print(' --> removing', dir_path)
                    shutil.rmtree(dir_path)
                else:
                    print(' --> skipping', dir_path)
            else:
                print(' --> not a halo output directory')

    def catalog_set_ext(self):
        known_types = ''
        for typ in self.catalog_known_ext:
            if typ == self.catalog_known_ext[- 1]:
                segment = typ + '.'
            else:
                segment = typ + ', '

            known_types += segment
        ext = tkSimpleDialog.askstring(
            title='Catalog Extension Type',
            prompt='Enter new Ext \nchoose from: ' + known_types,
            parent=self.plotwindow,
            initialvalue=self.catalog_extention)
        self.catalog_extention = ext
    def catalog_set_STIPS(self):
        pass

    def catalog_set_fits(self):
        pass

    def catalog_set_hdf5(self):
        pass

    def catalog_set_txt(self):
        pass

    def catalog_set_csv(self):
        pass

    def catalog_makeall(self):

        # If there are no regions to use.
        if not self.regions.keys():
            print('there are no regions to catalog')
            return

        # Make a catalog for each halo (halo) in self.regions.
        for halo in self.regions.keys():

            # Make a record table for this halos regions (rec_table).
            rec_table = self.record_table()

            # step through this halos regions (region).
            for region in self.regions[halo][2:]:

                # Terminal.
                print('loading region:', region['name'], 'for', halo)
                for key in region.keys():
                    print(' -> ', key, ' : ', region[key])

                row = output_catalog(
                    region, self.catalog_type, self.catalog_dir)
                rec_table.add_row(row)
                print(region['name'], 'done\n')

            print('writing record log text file')
            for halo in rec_table['halo']:
                table_fh = os.path.join(
                    self.catalog_dir, halo, 'recordtable.txt')
                with open(table_fh, 'w') as record_text:
                    for line in rec_table.pformat(max_width=200):
                        record_text.write(line + '\n')

            print('saving reference plot to local catalog file')
            for plot_fh in os.listdir(self.plot_dir):
                if plot_fh[:6] == halo:
                    src = os.path.join(self.plot_dir, plot_fh)
                    dest = os.path.join(self.catalog_dir, halo, plot_fh)
                    print('moving', plot_fh)
                    print('from:', src)
                    print('to:', dest)
                    shutil.move(src, dest)

    def make_cmd(self):
        gplt.plot_cmd(self)

    def plot_set_show_colorbar(self):
        if self.show_cb:
            self.show_cb = False
        else:
            self.show_cb = True

    def plot_show_colorbar(self):
        if not self.plot_has_cb:
            print('no color bar to plot')
            return
        image = PImage.open(self.plot_cb_fh)
        photo = PhotoImage(image)
        self.plot_cb = Label(master=self.plotwindow, image=photo)
        self.plot_cb.image = photo
        self.plot_cb.grid(row=3, column=0)

    def plot_set_ext_png(self):
        self.plot_extention = '.png'

    def plot_set_ext_pdf(self):
        self.plot_extention = '.pdf'

    def plot_set_ext_jpeg(self):
        self.plot_extention = '.jpeg'

    def plot_set_dpi_100(self):
        self.plot_dpi = 100

    def plot_set_dpi_200(self):
        self.plot_dpi = 200

    def plot_set_dpi_400(self):
        self.plot_dpi = 400

    def plot_set_dpi_600(self):
        self.plot_dpi = 600

    def plot_set_dpi_800(self):
        self.plot_dpi = 800

    def plot_set_radius_150(self):
        self.plot_radius_kpc = 150

    def plot_set_radius_200(self):
        self.plot_radius_kpc = 200

    def plot_set_radius_250(self):
        self.plot_radius_kpc = 250

    def plot_set_radius_300(self):
        self.plot_radius_kpc = 300

    def plot_set_directory(self):
        new_dir = tkFileDialog.askdirectory(
            title="Select New Plot Directory",
            mustexist=True,
            parent=self.master,
            initialdir=os.path.curdir)
        if os.path.isdir(new_dir):
            self.plot_dir = new_dir

    def plot_set_units_kpc(self):
        self.plot_units = 'kpc'

    def plot_set_units_arcsec(self):
        self.plot_units = 'arcsec'

    def plot_set_units_arcmin(self):
        self.plot_units = 'arcmin'

    def plot_set_units_deg(self):
        self.plot_units = 'degree'

    def plot_remove_all(self):
        for plot_file in os.listdir(self.plot_dir):
            if plot_file.endswith(self.plot_extention):
                os.remove(os.path.join(self.plot_dir, plot_file))
                print(' --> removed', plot_file)

    def plot_halo(self):
        if self.has_plot:
            self.clear_label_image()
        else:
            self.plot = Label(master=self.plotwindow, image='', cursor='cross')
        print('making plot file handle')
        self.plot_fh = os.path.join(self.plot_dir, self.halo + '.png')
        print(self.plot_fh)
        print('done')
        print('plotting')
        gplt.plot_halo(self)
        print('opening plot:', self.plot_fh)
        image = PImage.open(self.plot_fh).resize(self.plot_size)
        print('done')
        print('converting to TKimage')
        photo = PhotoImage(image)
        print('done')
        print('making label for plot display')
        self.plot = Label(master=self.plotwindow, image=photo, cursor='cross')
        print('done')
        print('attaching image to label')
        self.plot.image = photo
        print('binding click function')
        self.plot.bind("<Button-3>", self.right_click)
        self.plot.bind("<Button-1>", self.mouse_location)
        print('packing image label')
        self.plot.grid(row=0, column=0)
        if self.plot_has_cb and self.show_cb:
            self.plot_show_colorbar()



    def h4rg_10(self):
        self.region_pixel_size = 0.11
        self.region_ccd_size = (4088, 4088)

    def region_set_ccd_size(self):
        ccd_x = tkSimpleDialog.askinteger(
            title='Set WFI Active Detector CCD Size X-Axis',
            prompt='Enter new X length in pixels\nNumber of active detector columns: ',
            minvalue=100,
            maxvalue=16000,
            parent=self.plotwindow,
            initialvalue=self.region_ccd_size[0])
        ccd_y = tkSimpleDialog.askinteger(
            title='Set WFI Active Detector CCD Size Y-Axis',
            prompt='Enter new Y length in pixels\nNumber of active detector rows ',
            minvalue=100,
            maxvalue=16000,
            parent=self.plotwindow,
            initialvalue=self.region_ccd_size[1])
        self.region_ccd_size = (ccd_x, ccd_y)

    def region_set_pix_size(self):
        pixel = tkSimpleDialog.askfloat(
            title='Set WFI Pixel Scale',
            prompt='Enter new WFI detector pixel scale in arcseconds/pixel ',
            minvalue=0.05,
            maxvalue=0.17,
            parent=self.plotwindow,
            initialvalue=self.region_pixel_size)
        self.region_pixel_size = pixel

    def region_remove(self):
        answer = tkSimpleDialog.askinteger(
            title='Remove Region',
            prompt='Which region do you want to delete?\nEnter the integer number as seen on the plot')
        del self.regions[str(self.halo)][int(answer) + 1]
        self.regions[str(self.halo)][1] -= 1
        self.plot_halo()
        self.update_status_bar()

    def place_region(self):
        gplt.make_box(self)
        gplt.plot_halo(self)
        self.plot_halo()
        self.update_status_bar()

    def right_click(self, event):
        self.mousexy_to_plotxy(event.x, event.y)
        self.rightclick.tk_popup(event.x_root, event.y_root, 0)

    def mousexy_to_plotxy(self, mouse_x, mouse_y):
        offset_x = 400
        offset_y = 403
        visable_area = self.plot_radius_kpc * 2.0
        scale_x = visable_area / 621
        scale_y = visable_area / 614
        bump_x = -5
        bump_y = -1
        self.current_x = round((((mouse_x - offset_x) * scale_x) + bump_x))
        self.current_y = round(
            (-1 * (((mouse_y - offset_y) * scale_y) + bump_y)))
        self.smsg_mouse_xy = '[mouse xy - x: ' + \
            str(mouse_x) + '  y: ' + str(mouse_y) + ']   '
        self.smsg_kpc_xy = '[KPC - x: ' + \
            str(self.current_x) + '  y: ' + str(self.current_y) + ']    '
        self.update_status_bar()

    def mouse_location(self, event):
        self.mousexy_to_plotxy(event.x, event.y)
        print(self.smsg_status)
        print(self.master.geometry())

    def make_plotwindow(self):
        self.plotwindow = PanedWindow(
            master=self.master,
            width=self.plotwindow_size,
            height=self.plotwindow_size,
            bd=10)
        self.plotwindow.grid(row=1, column=0)
        self.plotwindow.bind("<Button-1>", self.mouse_location)

    def clear_label_image(self):
        self.plot.destroy()
        if self.show_cb and self.plot_has_cb:
            try:
                self.plot_cb.destroy()
            except AttributeError as e:
                print(e, file=os.sys.stderr)

    def update_status_bar(self):
        self.smsg_region = '[current regions: ' + \
            str(self.regions[str(self.halo)][1]) + ']    '
        self.smsg_status = self.smsg_region + self.smsg_mouse_xy + self.smsg_kpc_xy
        self.statusbar.set(self.smsg_status)

    def load_halo(self):
        self.load_halo_grid()
        self.load_halo_table()
        self.plot_halo()
        self.has_plot = True
        if self.plot_has_cb and self.show_cb:
            self.plot_show_colorbar()

    def load_halo_table(self):
        for fh in os.listdir(self.table_dir):
            comps = fh.split('_')
            if not comps[0] == self.halo:
                continue
            if fh.split('_')[2] == self.filter_type:
                self.table_fh = os.path.join(self.table_dir, fh)
                self.mag_table_fh = os.path.join(
                    self.table_dir, fh.replace('h158', 'CMD'))
                self.pxpy_table_fh = os.path.join(
                    self.table_dir, fh.replace('h158', 'pxpy'))
                print(self.table_fh)
                print(self.mag_table_fh)
                print(self.pxpy_table_fh)
                break
        self.table = load_table(self.table_fh)
        self.distance = self.table.meta['d_mpc']
        self.filter = self.table.meta['f_type']
        self.n_sats = self.table.meta['n_sats']
        self.sat0 = self.table.meta['sat0']
        self.sat1 = self.table.meta['sat1']
        self.abs_limit = self.table.meta['abm_lim']
        self.table.pprint()

    def load_halo_grid(self):
        self.grid_fh = tkFileDialog.askopenfilename(
            initialdir=self.grid_dir,
            title="Select new halo grid file")
        self.halo = str(os.path.split(self.grid_fh)[-1].split('_')[0])
        if not str(self.halo) in self.regions.keys():
            print('starting new region dict for', self.halo)
            self.regions[str(self.halo)] = [0, 0]

        self.grid = np.load(self.grid_fh)
        self.grid_center_x = self.grid.shape[1] / 2.0
        self.grid_center_y = self.grid.shape[0] / 2.0
        print(self.halo)
        self.master.title('StarCat - ' + self.halo)


def run_gui():
    root = Tk()
    app = App(root)
    root.mainloop()

if __name__ == '__main__':
    run_gui()
