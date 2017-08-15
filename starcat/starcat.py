from __future__ import division, absolute_import, print_function

import os

import printlib as plib

all_menus = {
    '1': 'menue name',
    '2': 'menue name',
    '3': 'menue name',
    '4': 'menue name'
}


def main_menu():

	# clear terminal
	if os.name == 'nt':
		clr_cmd = 'cls'
	else:
		clr_cmd = 'clear'
	os.system(clr_cmd)

	print(plib.main_title)
	print(plib.main_menu_title)

	for key in all_menus.keys():

		print('[', key, ']  --> ',all_menus[key])


