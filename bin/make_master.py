#!/usr/bin/env python

import glob
import numpy as np
import os
import sys
import time
import tqdm

from astropy.io import fits
from multiprocessing import Pool

from lyacolore import catalog, parse, utils

################################################################################

#Script to produce a master file from CoLoRe's output files.

################################################################################

"""
Set up the file locations and filename structures.
Also define option preferences.
"""

args = parse.get_args(sys.argv)

################################################################################

#Define the input file names and numbers.
input_files = utils.get_in_file_names(args.in_dir,args.input_filename_prefix)
if len(input_files) == 0:
    raise ValueError('No files found at {}'.format(args.in_dir))

#Get the simulation parameters from the parameter file.
simulation_parameters = utils.get_simulation_parameters(args.in_dir,args.param_file)

################################################################################
"""
Produce the data required to make a master file using multiprocessing (one task per original file).
Join the outputs from the many processes together.
Save the master file, and a similarly structured file containing QSOs with 'bad coordinates'.
"""

print('Working on master data...')
start = time.time()

master_data, bad_coordinates_data, cosmology_data = catalog.make_master_data(input_files,args,apply_footprint=True)
catalog.save_master_data(master_data,bad_coordinates_data,cosmology_data,args)

# TODO: this should really go in make_transmission. Issue is with multi-node bash script - it needs the directories to work out which pixels we have and thus how to distribute the work.

print('\nCreating the output file structure...')
#Make the new file structure
pixel_list = list(sorted(set(master_data['PIXNUM'])))
utils.make_file_structure(args.out_dir,pixel_list)
print(' -> Done!')
