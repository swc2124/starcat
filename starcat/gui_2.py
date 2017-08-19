'''
===========================================================================
starcat v0.1.0
September 2016
Sol W. Courtney
swc2124@columbia.edu
Columbia University NYC, NY
===========================================================================
this is the main file to be run
python starcat.py
'''
from __future__ import division, absolute_import, print_function
import os
import shutil
import sys
import pickle
import numpy as np
from PIL import Image as PImage, ImageTk
from Tkinter import *
import tkFileDialog
import tkSimpleDialog

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from astropy.table import Table
import astropy.units as u
import gui_plotlib as gplt

# tkSimpleDialog.askinteger(title, prompt [,options])
# tkSimpleDialog.askfloat(title, prompt [,options])
# tkSimpleDialog.askstring(title, prompt [,options])


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

        # =====================================================================
        #                              SETUP
        # =====================================================================

        # main window
        self.master = master
        self.master.title('StarCat')
        self.master.resizable(True, True)

        # file tree
        self.data_dir = 'not selected'
        for pth in os.path.abspath(os.path.curdir).split(os.path.sep)[:-1]:
            if self.data_dir == 'not selected':
                self.data_dir = os.path.join(pth, os.path.sep)
            self.data_dir = os.path.join(self.data_dir, pth)
        self.data_dir = os.path.join(self.data_dir, 'data')
        self.save_dir = os.path.join(self.data_dir, 'saved')
        self.table_dir = os.path.join(self.data_dir, 'tables')
        self.plot_dir = os.path.join(self.data_dir, 'plots')
        self.catalog_dir = os.path.join(self.data_dir, 'catalogs')
        self.grid_dir = os.path.join(self.data_dir, 'grids')

        # regions
        self.regions = []
        self.number_of_regions = 0
        self.number_of_regions_tot = 0
        self.region_size = 10

        # status messages
        self.smsg_region = '[current regions: ' + \
            str(self.number_of_regions) + ']    '
        self.smsg_mouse_xy = '[] '
        self.smsg_kpc_xy = '[]   '

        # plot display
        self.plot_size = (800, 800)
        self.plotwindow_size = 1000
        self.has_plot = False
        master.minsize(width=self.plot_size[
                       0] + 50, height=self.plot_size[0] + 100)

        # plot settings
        self.plot_dpi = 400
        self.plot_units = 'kpc'
        self.plot_known_units = ['kpc', 'arcmin', 'arcsec', 'degree']
        self.plot_extention = '.png'
        self.plot_radius_kpc = 150

        # catalog settings
        self.catalog_extention = '.fits'

        # =====================================================================
        #                              WINDOWS
        # =====================================================================

        # status bar
        self.statusbar = StatusBar(self.master)
        self.statusbar.grid(row=2, column=0, columnspan=2)

        # plot window
        self.make_plotwindow()

        # =====================================================================
        #                               MENUS
        # =====================================================================

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
            label="Set plot dpi",
            command=self.plot_set_dpi)

        plot_menu.add_command(
            label="Set plot directory",
            command=self.plot_set_directory)

        plot_menu.add_command(
            label="Set plot units",
            command=self.plot_set_units)
        plot_menu.add_command(
            label="Set plot radius",
            command=self.plot_set_radius)

        plot_menu.add_command(
            label="Set plot extension",
            command=self.plot_set_ext)
        plot_menu.add_separator()
        plot_menu.add_command(
            label="Remove all",
            command=self.plot_remove_all,
            background='red')

        topmenu.add_cascade(
            label="Plot",
            menu=plot_menu)

        # top menu - catalog menu
        self.catalog_menu = Menu(topmenu)
        self.catalog_menu.add_command(
            label='Make Catalogs',
            command=self.catalog_makeall)
        self.catalog_menu.add_separator()
        self.catalog_menu.add_command(
            label='Remove all',
            command=self.catalog_remove_all,
            background='red')
        topmenu.add_cascade(
            label='Catalog',
            menu=self.catalog_menu)

        # right click menu
        self.rightclick = Menu(self.plotwindow, tearoff=0)
        self.rightclick.add_command(
            label='Select Region',
            command=self.place_region)
        self.rightclick.add_command(
            label='Remove region',
            command=self.region_remove)

    # =====================================================================
    #                              FUNCTIONS
    # =====================================================================
    def plot_set_ext(self):
        print('[plot_set_ext]')
        print('current ext:', self.plot_extention)
        ext = tkSimpleDialog.askstring(
            title='Plot Ext',
            prompt='Enter new Extention',
            parent=self.plotwindow,
            initialvalue=self.plot_extention)
        self.plot_extention = ext
        print('new ext:', self.plot_extention)
        print('[plot_set_ext complete]')

    def plot_set_dpi(self):
        print('[plot_set_dpi]')
        print('current dpi:', self.plot_dpi)
        dpi = tkSimpleDialog.askinteger(
            title='Plot Dpi',
            prompt='Enter new dpi',
            minvalue=50,
            maxvalue=1000,
            parent=self.plotwindow,
            initialvalue=self.plot_dpi)
        self.plot_dpi = dpi
        print('new dpi:', self.plot_dpi)
        print('[plot_set_dpi complete]')

    def plot_set_radius(self):
        print('[plot_set_radius]')
        rad = tkSimpleDialog.askinteger(
            title='Plot Radius Kpc',
            prompt='Enter new radius in Kpc',
            minvalue=50,
            maxvalue=300,
            parent=self.plotwindow,
            initialvalue=self.plot_radius_kpc)
        self.plot_radius_kpc = rad
        print('[plot_set_radius complete]')

    def plot_set_directory(self):
        print('[plot_set_directory]')
        print('current plot directory:', self.plot_dir)
        new_dir = tkFileDialog.askdirectory(
            title="Select New Plot Directory",
            mustexist=1)
        if os.path.isdir(new_dir):
            self.plot_dir = new_dir
        print('new plot directory:', self.plot_dir)
        print('[plot_set_directory complete]')

    def plot_set_units(self):
        print('[plot_set_units]')
        print('current units', self.plot_units)
        known_types = ''
        for typ in self.plot_known_units:
            if typ == self.plot_known_units[- 1]:
                segment = typ + '.'
            else:
                segment = typ + ', '

            known_types += segment
        unit = tkSimpleDialog.askstring(
            title='Plot Units',
            prompt='Enter units \nchoose from: ' + known_types,
            parent=self.plotwindow,
            initialvalue=self.plot_units)
        if unit not in self.plot_known_units:
            print('unknown unit', unit)
            return
        self.plot_units = unit
        print('new units', self.plot_units)
        print('[plot_set_units complete]')

    def plot_remove_all(self):
        print('[plot_remove_all]')
        for plot_file in os.listdir(self.plot_dir):
            if plot_file.endswith(self.plot_extention):
                os.remove(os.path.join(self.plot_dir, plot_file))
                print(' --> removed', plot_file)
        print('[plot_remove_all complete]')

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

    def catalog_remove_all(self):
        print('[catalog_remove_all]')
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
        print('[catalog_remove_all complete]')

    def catalog_makeall(self):
        print('[catalog_makeall]')

        if not self.regions:
            print('there are not regions to catalog')
            return
        record_table = self.record_table()
        for region in self.regions:
            print('loading region:', region['name'])
            for key in region.keys():
                print(' -> ', key, ' : ', region[key])
            catalog_name = region['name']
            catalog_dir = os.path.join(self.catalog_dir, region['halo'])

            # TODO - ask to overwrite existing data
            # ===============================================================================
            # TODO
            if not os.path.isdir(catalog_dir):
                print('making catalog directory')
                os.mkdir(catalog_dir)
                print(catalog_dir)
            else:
                print('catalog directory exists')
                print('contents:')
                for f in os.listdir(catalog_dir):
                    print('  --> ', f)

            print('making catalog file handle')
            catalog_fh = os.path.join(
                catalog_dir, region['name'] + self.catalog_extention)
            print(catalog_fh)
            if os.path.isfile(catalog_fh):
                print('catalog exists')
                try:
                    number = '_0' + \
                        str(int(catalog_fh.split('.')[-2][-2:]) + 1)
                except:
                    number = '_01'
                # TODO - ask to overwrite existing data
                # ===============================================================================
                # TODO
                catalog_fh = os.path.join(
                    catalog_dir, region['name'] + number + self.catalog_extention)
            print('finding stars within region boundaries')
            print('x0:', region['x0'])
            print('x1:', region['x1'])
            print('y0:', region['y0'])
            print('y1:', region['y1'])
            idx = np.nonzero(
                np.logical_and(
                    np.logical_and(
                        self.table['x_int'] >= region['x0'],
                        self.table['x_int'] <= region['x1']),
                    np.logical_and(
                        self.table['y_int'] >= region['y0'],
                        self.table['y_int'] <= region['y1'])
                )
            )[0]
            print('found', len(idx), 'stars within region')
            print('making catalog from table')
            fits_table = self.table[idx]
            fits_table.pprint()
            # fits_table.meta['']
            if 'spinbin_output_fh' in fits_table.meta.keys():
                del fits_table.meta['spinbin_output_fh']
            print('writing catalog to disc')
            fits_table.write(catalog_fh, format='fits')
            print('adding row to record table')
            row = [region['halo'], region['name'], len(idx), fits_table['mact'].sum(),
                   region['x0'], region['x1'], region['y0'], region['y1'],
                   fits_table['feh'].min(), fits_table[
                'feh'].mean(), fits_table['feh'].max(),
                fits_table['age'].min(), fits_table['age'].mean(), fits_table['age'].max()]
            record_table.add_row(row)
            print(region['name'], 'done\n')

        print('writing record log text file')
        for halo in record_table['halo']:
            table_fh = os.path.join(self.catalog_dir, halo, 'recordtable.txt')
            with open(table_fh, 'w') as record_text:
                for line in record_table.pformat(max_width=200):
                    record_text.write(line + '\n')
        print('[catalog_makeall complete]')

    def region_remove(self):
        answer = tkSimpleDialog.askinteger(
            title='Remove Region',
            prompt='Which region do you want to delete?\nEnter the integer number as seen on the plot')
        del self.regions[int(answer) - 1]
        self.number_of_regions -= 1
        self.plot_halo()
        self.update_status_bar()
        print('remove region', answer)

    def regions_display(self):
        region_list = ['region ' + str(i) for i in range(0, len(self.regions))]

    def make_plotwindow(self):
        print('[make_plotwindow]')
        self.plotwindow = PanedWindow(
            master=self.master,
            width=self.plotwindow_size,
            height=self.plotwindow_size)
        self.plotwindow.grid(row=1, column=0, columnspan=2)
        self.plotwindow.bind("<Button-1>", self.mouse_location)
        print('[make_plotwindow complete]')

    def clear_label_image(self):
        print('[clear_label_image]')
        print('self.has_plot:', self.has_plot)
        # self.plotwindow.forget(self.plot)
        # self.plot.config(image='')
        # self.master.forget(self.plot)
        # self.make_plotwindow()
        self.plot.destroy()
        print('[clear_label_image complete]')

    def update_status_bar(self):
        print('[update_status_bar]')
        self.smsg_region = '[current regions: ' + \
            str(self.number_of_regions) + ']    '
        self.smsg_status = self.smsg_region + self.smsg_mouse_xy + self.smsg_kpc_xy
        self.statusbar.set(self.smsg_status)
        print('[update_status_bar complete]')

    def place_region(self):
        print('[place_region]')
        gplt.make_box(self)
        gplt.plot_halo(self)
        self.plot_halo()
        self.update_status_bar()
        print('[place_region complete]')

    def right_click(self, event):
        # saves mouse xy in plot numbers
        # sends to rightclick menu
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

    def load_halo(self):
        print('[load_halo]')
        self.load_halo_grid()
        self.load_halo_table()
        self.plot_halo()
        self.has_plot = True
        print('[load_halo complete]')

    # def display_table(self):
    #    table = Label(master=self.datawindow, text=self.table.pformat(max_lines=100), anchor=W)
    #    table.grid()

    def load_halo_table(self):
        print('[load_halo_table]')
        for fh in os.listdir(self.table_dir):
            if fh.startswith(self.halo):
                self.table_fh = os.path.join(self.table_dir, fh)
                print(self.table_fh)
                break
        self.table = Table.read(self.table_fh, format='hdf5', path='data')
        self.distance = self.table.meta['d_mpc']
        self.filter = self.table.meta['f_type']
        self.n_sats = self.table.meta['n_sats']
        self.sat0 = self.table.meta['sat0']
        self.sat1 = self.table.meta['sat1']
        self.abs_limit = self.table.meta['abm_lim']
        self.table.pprint()
        print('[load_halo_table complete]')

    def load_halo_grid(self):
        print('[load_halo_grid]')
        self.grid_fh = tkFileDialog.askopenfilename(
            initialdir=self.grid_dir,
            title="Select new halo grid file")
        self.halo = os.path.split(self.grid_fh)[-1].split('_')[0]
        self.grid = np.load(self.grid_fh)
        print(self.halo)
        # self.current_halo_bar.set(str(self.halo))
        self.master.title('StarCat - ' + self.halo)
        # self.plot_halo()
        print('[load_halo_grid complete]')
    '''
    def plot_halo(self):
        if self.dataplot:
            del self.dataplot
        self.dataplot = FigureCanvasTkAgg(gplt.plot_halo(self), master=self.plotwindow)
        self.dataplot.show()
        self.dataplot.get_tk_widget().grid()

    '''

    def plot_halo(self):
        print('[plot_halo]')
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
        '''
        if not self.halo + '.png' in os.listdir(self.plot_dir):
            print('no plot exists')
            print('plot_fh: ', self.plot_fh)
            print('plot_dir :',)
            for f in os.listdir(self.plot_dir):
                print('  --> ', f)
            gplt.plot_halo(self)
        else:
            print('loading existing plot')
        '''
        print('opening plot:', self.plot_fh)
        image = PImage.open(self.plot_fh).resize(self.plot_size)
        print('done')
        print('converting to TKimage')
        photo = ImageTk.PhotoImage(image)
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
        print('[plot_halo complete]')


root = Tk()
#top = Toplevel()

app = App(root)
root.mainloop()
root.destroy()
