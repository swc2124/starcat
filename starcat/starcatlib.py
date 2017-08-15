from __future__ import division, absolute_import, print_function


import os
import numpy as np
from astropy.table import Table
from astropy.table import Column
from astropy import units as u
from astropy.units.function.units import dex
import ebf
from c_functions import integerize
from matplotlib import pyplot as plt


def make_boxes(_feh, _age, spread=0.025, feh_spreadbump=0.0, age_spreadbump=0.0):
    '''
    box coord set making function
    make a list of 4 boxes from feh min/max and age min/max
    boxlist = make_boxes()
    returns [[(3.85, 4.15), (-1.25, -0.95)],
             [(3.85, 4.15), (-0.85, -0.55)],
             [(10.85, 11.15), (-1.25, -0.95)],
             [(10.85, 11.15), (-0.85, -0.55)]]
    '''
    # Fe/H regions : fehrange0 = (f_start0, f_stop0), fehrange1 = (f_start1,
    # f_stop1)
    feh_spread = ((_feh[1] - _feh[0]) * spread) + feh_spreadbump
    feh_0 = _feh[0] - feh_spread, _feh[0] + feh_spread
    feh_1 = _feh[1] - feh_spread, _feh[1] + feh_spread
    feh_regions = [feh_0, feh_1]

    # age regions : agerange0 = (a_start0, a_stop0), agerange1 = (a_start1,
    # a_stop1)
    age_spread = ((_age[1] - _age[0]) * spread) + age_spreadbump
    age_0 = _age[0] - age_spread, _age[0] + age_spread
    age_1 = _age[1] - age_spread, _age[1] + age_spread
    age_regions = [age_0, age_1]

    # make list of boxes for place_box() function
    box_list = []
    for _Age in age_regions:
        for _Feh in feh_regions:
            box_list.append([(_Age[0], _Age[1]), (_Feh[0], _Feh[1])])

    # return list
    return box_list


def place_box(_ax, _box, r_number):
    '''
    box placing function
    make a single box from a set of numbers like this: box = ((x0, x1), (y0, y1))
    _ax.vlines(x0, y0, y1)
    _ax.vlines(x1, y0, y1)
    _ax.hlines(y0, x0, x1)
    _ax.hlines(y1, x0, x1)
    '''
    # unpack box tuples
    _age, _feh = _box

    # find box center
    a_center = _age[0] + ((_age[1] - _age[0]) * 0.4)
    f_center = _feh[1] + ((_feh[1] - _feh[0]) * 0.7)

    # set box label
    box_label = 'R ' + str(r_number)
    _ax.text(a_center, f_center, box_label,
             ha='center',
             fontsize=25,
             color='k',
             bbox={'facecolor': 'white', 'alpha': 0.75, 'pad': 2})

    # params
    lwidth = 3
    lstyle = 'dashed'
    lclr = 'r'

    # set lines
    _ax.vlines(_age[0], _feh[0], _feh[1],
               colors=lclr,
               linestyles=lstyle,
               linewidth=lwidth)
    _ax.vlines(_age[1], _feh[0], _feh[1],
               colors=lclr,
               linestyles=lstyle,
               linewidth=lwidth)
    _ax.hlines(_feh[0], _age[0], _age[1],
               colors=lclr,
               linestyles=lstyle,
               linewidth=lwidth)
    _ax.hlines(_feh[1], _age[0], _age[1],
               colors=lclr,
               linestyles=lstyle,
               linewidth=lwidth)


def make_region_plot(_feh, _age,
                     _ebf_feh, _ebf_age,
                     feh_range,
                     age_range,
                     d_cut=100,
                     _halo=halo):
    '''
    age vs Metallicity plot with selected regions in red boxes
    EBF and DAT stellar data plotted over each other.
    Make and save a plot of age vs feh with the selected regions boxed in red.
    _feh, _age: array or list
    _ebf_feh, _ebf_age : EBF data --> ebf_data['feh'], ebf_data['age']
    feh_range, age_range : tupple --> (feh_0, feh_1), (age_0, age_1) center of regions.
    '''
    # figure & figsize
    fig = plt.figure(figsize=(10, 10))

    # subplot
    ax = fig.add_subplot(111)

    # plot EBF and Dat data over each other
    ax.scatter(np.power(10, _ebf_age[::d_cut]) / 1e9,
               _ebf_feh[::d_cut],
               s=8,
               c='k',
               alpha=.15)

    ax.scatter(age[::d_cut],
               feh[::d_cut],
               s=4,
               c='grey',
               alpha=.15)

    # plot boxes one at a time
    for i, box in enumerate(make_boxes(_feh=feh_range, _age=age_range)):
        place_box(ax, box, (i + 1))

    # plot limits
    ax.set_xlim([1.0, 12.7])
    ax.set_ylim([-1.5, -0.4])

    # plot labels
    ax.set_title(
        'EBF (black) & Dat (grey) Data\nSelected Ages and Metallicities\n' + _halo)
    ax.set_xlabel('Age (Gyr)')
    ax.set_ylabel('Fe/H (Dex)')

    # plot grid
    plt.grid()

    # save plot or show plot
    # plt.show()
    fig.savefig('./' + _halo + '_ebf_and_dat_regions')
    plt.close()


def make_mt_table():
    '''
    helper function
    returns an empty astropy table
    Returns:
        [astropy.table] -- [description]
    '''
    Names = ['region', 'feh_min', 'feh_mean', 'feh_max',
             'log_mass_tot', 'mass_min', 'mass_mean', 'mass_max',
             'age_min', 'age_mean', 'age_max']
    Dtype = ['i', 'f', 'f', 'f',
             'f', 'f', 'f', 'f',
             'f', 'f', 'f']
    return Table(names=Names, dtype=Dtype)


def make_dat_table(age, feh, feh_range, age_range):
    '''make a table showing the stats for .dat file data

    [description]

    Arguments:
        age {[array]} -- [description]
        feh {[array]} -- [description]
        feh_range {[tupple]} -- [description]
        age_range {[tupple]} -- [description]

    Returns:
        [table] -- [description]
    '''
    # make an empty table
    dat_table = make_mt_table()

    # itterate through the boxed regions
    for i, age_feh in enumerate(make_boxes(feh_range, age_range)):

        # unpack limits for this box
        _age, _feh = age_feh
        age0, age1 = _age
        feh0, feh1 = _feh

        # get indices for this box
        idx = np.nonzero(
            np.logical_and(
                np.logical_and(
                    age >= age0,
                    age < age1),
                np.logical_and(
                    feh >= feh0,
                    feh < feh1)
            ))[0]

        # skip if empty (should not happen)
        if not len(idx):
            print idx
            continue

        # add this box's stats to table
        dat_table.add_row([
            i + 1,
            round(feh[idx].min(), 2),
            round(feh[idx].mean(), 2),
            round(feh[idx].max(), 2),
            round(np.log10(mass[idx].sum()), 2),
            round(mass[idx].min(), 2),
            round(mass[idx].mean(), 2),
            round(mass[idx].max(), 2),
            round(age[idx].min(), 2),
            round(age[idx].mean(), 2),
            round(age[idx].max(), 2)])

    print('-----------------')
    print('Dat file results')
    print('-----------------')
    dat_table.pprint(max_width=200)
    return dat_table


def make_ebf_table_and_cmd(ebf_data, feh_range, age_range, _halo=halo):

    dat_table = make_dat_table(age, feh, feh_range, age_range)

    # new table for EBF data
    ebf_table = make_mt_table()
    ebf_table.add_column(Column(name='N_stars'))
    # new figure for EBF CMD plot
    fig = plt.figure(figsize=(30, 30))

    # iterate through the regions
    for i, age_feh in enumerate(make_boxes(feh_range, age_range)):

        # unpack limits
        age0, age1 = age_feh[0]
        feh0, feh1 = age_feh[1]

        # get indices
        idx = np.nonzero(
            np.logical_and(
                np.logical_and(
                    np.power(10, ebf_data['age']) / 1e9 >= age0,
                    np.power(10, ebf_data['age']) / 1e9 < age1),
                np.logical_and(
                    ebf_data['feh'] >= feh0,
                    ebf_data['feh'] < feh1)
            ))[0]

        # color array for EBF CMD
        color = ebf_data['wfirst-hst_z087'][idx] - \
            ebf_data['wfirst-hst_h158'][idx]

        # new subplot for this region
        ax = fig.add_subplot(111)

        # Legend label
        L_0 = 'R ' + str(i + 1)
        L_1 = ' | ' + str(round(feh0, 2)) + '/' + str(round(feh1, 2)) + ' Dex'
        L_2 = ' - ' + str(round(age0, 2)) + '/' + str(round(age1, 2)) + ' Gyr'
        L_3 = ' | ' + str(len(idx)) + ' stars'
        Label = L_0 + L_1 + L_2 + L_3

        # add this region's stars to the EBF CMD plot
        ax.scatter(color,
                   ebf_data['wfirst-hst_z087'][idx],
                   s=15,
                   alpha=1,
                   label=Label)

        # add row to EBF table
        ebf_table.add_row([i + 1,
                           round(ebf_data['feh'][idx].min(), 2),
                           round(ebf_data['feh'][idx].mean(), 2),
                           round(ebf_data['feh'][idx].max(), 2),
                           round(np.log10(ebf_data['mact'][idx].sum()), 2),
                           round(ebf_data['mact'][idx].min(), 2),
                           round(ebf_data['mact'][idx].mean(), 2),
                           round(ebf_data['mact'][idx].max(), 2),
                           round(np.power(10, ebf_data[
                                 'age'][idx].min()) / 1e9, 2),
                           round(np.power(10, ebf_data['age'][
                                 idx].mean()) / 1e9, 2),
                           round(np.power(10, ebf_data[
                                 'age'][idx].max()) / 1e9, 2),
                           len(idx)])

    # x and y labels
    plt.title('Color Magnitude Diagram\n' + halo, fontsize=40)
    plt.xlabel('Z087 - H158', fontsize=40)
    plt.ylabel('Z087 AB_mag', fontsize=40)

    # x and y limits
    plt.xlim([-1.2, 3.2])
    plt.ylim([2.0, -4.5])

    # plot legend
    plt.legend(fontsize=30, markerscale=5, fancybox=True,
               shadow=True, framealpha=0.5)

    # plot tick sizes
    plt.tick_params(labelsize=30)

    # plot grid
    plt.grid()

    # save the plot
    fig.savefig('./' + _halo + '_cmd')

    # print some pretty output
    print '-----------------'
    print 'Dat file results'
    print '-----------------'
    dat_table.pprint(max_width=200)
    print '\n-----------------'
    print 'EBF file results'
    print '-----------------'
    ebf_table.pprint(max_width=200)

    return ebf_table


def make_fits_catalog(ebf_data, idx, table_name, halo=halo):

    table = Table()

    dtypes = [
        ('teff', u.Kelvin), ('wfirst-hst_f606w', u.mag),
        ('wfirst-hst_w149', u.mag), ('pz', u.kiloparsec),
        ('px', u.kiloparsec), ('py', u.kiloparsec),
        ('wfirst-hst_f160w', u.mag), ('feh', u.dex),
        ('wfirst-hst_y106', u.mag), ('wfirst-hst_f110w', u.mag),
        ('wfirst-hst_h158', u.mag), ('wfirst-hst_f184', u.mag),
        ('lum', u.solLum), ('mact', u.solMass),
        ('wfirst-hst_f814w', u.mag), ('wfirst-hst_f475w', u.mag),
        ('alpha', u.dex), ('wfirst-hst_z087', u.mag),
        ('wfirst-hst_f555w', u.mag), ('age', u.gigayear),
        ('wfirst-hst_j129', u.mag), ('smass', u.solMass)]

    for key, _unit in dtypes:

        _name = key
        if key.startswith('wfirst'):
            _name = key.replace('wfirst-hst_', '')

        new_col = Column(data=ebf_data[key][idx], name=_name, unit=_unit)
        table.add_column(new_col)
        # print '   --> ' + key

    # table.pprint()
    table.write('./' + table_name + '.fits', format='fits', overwrite=True)
    # print 'table saved'


# In[8]:

FEH_Range = (-1.215, -0.919)
AGE_Range = (3.0, 11.5)


# In[10]:

halos = [f.split('.')[0] for f in os.listdir('./') if f.endswith('.dat')]
halos.reverse()
halos


# In[11]:

for halo in halos:
    print halo
    age, feh = np.loadtxt('./' + halo + '.dat', unpack=True, usecols=(10, 11))
    ebf_data = ebf.read('C:\Users\\swc21\\Desktop\\Halos - Copy\\' + halo + '.ebf')
    print '  loaded'
    make_region_plot(feh, age, ebf_data['feh'], ebf_data['age'],
                     feh_range=FEH_Range,
                     age_range=AGE_Range,
                     _halo=halo)
    print '  done'


# In[ ]:

# make ebf table and CMD for ebf data
ebf_table = make_ebf_table_and_cmd(
    ebf_data, feh_range=FEH_Range, age_range=AGE_Range)

# set the EBF table column values for units
ebf_table['region'].unit = 'int'
ebf_table['feh_min'].unit = 'Dex'
ebf_table['feh_mean'].unit = 'Dex'
ebf_table['feh_max'].unit = 'Dex'
ebf_table['log_mass_tot'].unit = 'Log10(' + u.solMass.to_string() + ')'
ebf_table['mass_min'].unit = u.solMass
ebf_table['mass_mean'].unit = u.solMass
ebf_table['mass_max'].unit = u.solMass
ebf_table['age_min'].unit = u.gigayear
ebf_table['age_mean'].unit = u.gigayear
ebf_table['age_max'].unit = u.gigayear
ebf_table['N_stars'].unit = 'int'

# write a txt file with basic stats for each of the 4 regions
with open('./' + halo + '_Region_info.txt', 'w') as file:
    for line in ebf_table.pformat(max_lines=200, max_width=200, show_name=True, show_unit=True, show_dtype=True, align='^'):
        file.write(line + '\n')


# In[ ]:

# make a fits table for each
for i, age_feh in enumerate(make_boxes(_feh=FEH_Range, _age=AGE_Range)):

    # unpack limits
    age0, age1 = age_feh[0]
    feh0, feh1 = age_feh[1]

    # get indices
    idx = np.nonzero(
        np.logical_and(
            np.logical_and(
                np.power(10, ebf_data['age']) / 1e9 >= age0,
                np.power(10, ebf_data['age']) / 1e9 < age1),
            np.logical_and(
                ebf_data['feh'] >= feh0,
                ebf_data['feh'] < feh1)
        ))[0]

    # write the table as a fits file
    make_fits_catalog(ebf_data, idx, halo + '_Region_' + str(i + 1))


# In[ ]:


# In[ ]:


# In[ ]:


# In[ ]:


# In[ ]:


# In[ ]:


# In[ ]:


# In[ ]:

pt


# In[ ]:

table.write()
'''
=========================== ==== ===== ============= ==========
           Format           Read Write Auto-identify Deprecated
=========================== ==== ===== ============= ==========
                      ascii  Yes   Yes            No           
               ascii.aastex  Yes   Yes            No           
                ascii.basic  Yes   Yes            No           
     ascii.commented_header  Yes   Yes            No           
                  ascii.csv  Yes   Yes            No           
                 ascii.ecsv  Yes   Yes            No           
           ascii.fast_basic  Yes   Yes            No           
ascii.fast_commented_header  Yes   Yes            No           
             ascii.fast_csv  Yes   Yes            No           
       ascii.fast_no_header  Yes   Yes            No           
             ascii.fast_rdb  Yes   Yes            No           
             ascii.fast_tab  Yes   Yes            No           
          ascii.fixed_width  Yes   Yes            No           
ascii.fixed_width_no_header  Yes   Yes            No           
 ascii.fixed_width_two_line  Yes   Yes            No           
                 ascii.html  Yes   Yes           Yes           
                 ascii.ipac  Yes   Yes            No           
                ascii.latex  Yes   Yes           Yes           
            ascii.no_header  Yes   Yes            No           
                  ascii.rdb  Yes   Yes           Yes           
                  ascii.rst  Yes   Yes            No           
                  ascii.tab  Yes   Yes            No           
                       fits  Yes   Yes           Yes           
                       hdf5  Yes   Yes           Yes           
                   jsviewer   No   Yes            No           
                    votable  Yes   Yes           Yes           
                     aastex  Yes   Yes            No        Yes
                        csv  Yes   Yes           Yes        Yes
                       html  Yes   Yes            No        Yes
                       ipac  Yes   Yes            No        Yes
                      latex  Yes   Yes            No        Yes
                        rdb  Yes   Yes            No        Yes
=========================== ==== ===== ============= ==========


'''


# In[ ]:


# In[ ]:


# In[ ]:


# In[ ]:


# In[ ]:

idx = np.nonzero(np.logical_and(
    np.logical_and(px > -300.0, px < 300.0),
    np.logical_and(py > -300.0, py < 300.0)))
Ix, Iy = integerize(px[idx].astype(np.float64), py[idx].astype(np.float64))
len(Ix)


# In[ ]:

grid = np.zeros((600, 600), dtype=np.float64)
for i in range(len(Ix)):
    if Ix[i] > 599.0 or Iy[i] > 599.0:
        continue
    grid[Ix[i], Iy[i]] += mass[i]


# In[ ]:

from matplotlib import pyplot as plt
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)

ax.pcolormesh(np.log10(grid), cmap=plt.cm.plasma, vmin=0., vmax=6.5)
ax.axes.grid(linestyle='dashed', alpha=.5)
ticks = [str(i) for i in range(-300, 301, 100)]
fig.savefig('./heatmap')
plt.show()


# In[ ]:


age_bins = [(i, i + 1) for i in range(4, 12, 1)]
for age_start, age_stop in age_bins:

    Feh = feh[np.nonzero(np.logical_and(age >= age_start, age < age_stop))]
    ebf_Feh = halo02['feh'][np.nonzero(np.logical_and(
        halo02['age'] >= age_start, halo02['age'] < age_stop))]
    if not len(ebf_Feh):
        print 'bunk'
        continue
    print ''
    print 'age', age_start, age_stop
    print Feh.min()
    print Feh.mean()
    print Feh.max()
    print ''
    print ebf_Feh.min()
    print ebf_Feh.mean()
    print ebf_Feh.max()
    print ''


# In[ ]:

feh_bins = [(-3.5, -3.0), (-3.0, -2.5), (-2.5, -2.0),
            (-2.0, -1.5), (-1.5, -1.0), (-1.0, -0.5), (-0.5, 0.0)]
for feh_start, feh_stop in feh_bins:

    Age = age[np.nonzero(np.logical_and(feh >= feh_start, feh < feh_stop))]
    ebf_Age = halo02['age'][np.nonzero(np.logical_and(
        halo02['feh'] >= feh_start, halo02['feh'] < feh_stop))]
    if not len(ebf_Age):
        print 'bunk'
        continue
    print ''
    print 'feh', feh_start, feh_stop
    print ''
    print 'dat'
    print Age.min()
    print Age.mean()
    print Age.max()
    print ''
    print 'ebf'
    print round(np.power(10, ebf_Age.min()) / 1e9, 2)
    print round(np.power(10, ebf_Age.mean()) / 1e9, 2)
    print round(np.power(10, ebf_Age.max()) / 1e9, 2)
    print ''


# In[ ]:


# In[ ]:


# In[ ]:


# In[ ]:


# In[ ]:


# In[ ]:


# In[ ]:

halo02 = ebf.read('../Halos - Copy/halo08.ebf')
print '[age comp] \n\nDat file\n---------'
print round(age.min(), 2)
print round(age.mean(), 2)
print round(age.max(), 2)
print '\n\nEBF data\n---------'
print round(np.power(10, halo02['age'].min()) / 1e9, 2)
print round(np.power(10, halo02['age'].mean()) / 1e9, 2)
print round(np.power(10, halo02['age'].max()) / 1e9, 2)
print '\n[feh comp] \n\nDat file\n---------'
print round(feh.min(), 2)
print round(feh.mean(), 2)
print round(feh.max(), 2)
print '\n\nEBF data\n---------'
print round(halo02['feh'].min(), 2)
print round(halo02['feh'].mean(), 2)
print round(halo02['feh'].max(), 2)
