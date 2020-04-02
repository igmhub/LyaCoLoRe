#!/usr/bin/env python

import configargparse
import glob
import numpy as np
import os
import sys
import time

from astropy.io import fits
from multiprocessing import Pool

from lyacolore import catalog, parse, utils

################################################################################

#Script to produce a master file from CoLoRe's output files.

# TODO: Get rid of the need to specify file numbers?

################################################################################

"""
Set up the file locations and filename structures.
Also define option preferences.
"""

args = parse.get_args(sys.argv)

################################################################################

# TODO: print to confirm the arguments. e.g. "DLAs will be added"

#Check the value of N_side required is a power of 2.
if np.log2(args.nside)-int(np.log2(args.nside)) != 0:
    print('nside must be a power of 2!')
else:
    N_pix = 12*args.nside**2

#Define the original file structure
input_filename_structure = 'out_srcs_s1_{}.fits' #file_number
input_files = glob.glob(args.in_dir+'/'+input_filename_structure.format('*'))
file_numbers = utils.get_file_numbers(args.in_dir,input_filename_structure,input_files)
input_format = 'gaussian_colore'

#Get the simulation parameters from the parameter file.
simulation_parameters = utils.get_simulation_parameters(args.in_dir,args.param_file)

################################################################################

"""
Define the multiprocessing tracking functions
"""

#Define a progress-tracking function.
def log_result(retval):

    results.append(retval)
    N_complete = len(results)
    N_tasks = len(tasks)
    utils.progress_bar(N_complete,N_tasks,start_time)

    return

#Define an error-tracking function.
def log_error(retval):
    print('Error:',retval)
    return

################################################################################

"""
Produce the data required to make a master file using multiprocessing (one task per original file).
Join the outputs from the many processes together.
Save the master file, and a similarly structured file containing QSOs with 'bad coordinates'.
"""

print('Working on master data...')
start = time.time()

#Choose the QSO filtering we want.
QSO_filter = utils.make_QSO_filter(args.footprint,N_side=args.nside)

#Define the process to make the master data.
def make_master_data(file_name,file_number,input_format,N_side,minimum_z=args.min_cat_z):

    seed = args.seed*10**5 + file_number
    file_number, ID_data, cosmology, file_pixel_map_element, MOCKID_lookup_element = catalog.get_ID_data(file_name,file_number,input_format,args.nside,minimum_z=minimum_z,downsampling=args.downsampling,QSO_filter=QSO_filter,pixel_list=args.pixels,seed=seed)

    return [file_number, ID_data, cosmology, file_pixel_map_element, MOCKID_lookup_element]

#Set up the multiprocessing pool parameters and make a list of tasks.
tasks = [(input_files[i],file_number,input_format,args.nside) for i,file_number in enumerate(file_numbers)]

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = args.nproc)
    results = []
    start_time = time.time()

    for task in tasks:
        pool.apply_async(make_master_data,task,callback=log_result,error_callback=log_error)

    pool.close()
    pool.join()

print('\nSaving the master files...')

#Join the multiprocessing results into 'master' and 'bad_coordinates' arrays.
master_data, bad_coordinates_data, cosmology_data, file_pixel_map, MOCKID_lookup = catalog.join_ID_data(results,args.nside)

#Write master and bad coordinates files.
master_filename = args.out_dir + '/master.fits'
catalog.write_ID(master_filename,args.nside,master_data,cosmology_data,overwrite=args.overwrite)
print(' -> Master file contains {} objects.'.format(master_data.shape[0]))

if bad_coordinates_data.shape[0] > 0:
    bad_coordinates_filename = args.out_dir + '/bad_coordinates.fits'
    catalog.write_ID(bad_coordinates_filename,args.nside,bad_coordinates_data,cosmology_data,overwrite=args.overwrite)
    print(' -> "Bad coordinates" file contains {} objects.'.format(bad_coordinates_data.shape[0]))

#If desired, write the DRQ files for picca xcf to deal with.
if args.add_picca_drqs:
    for RSD_option in ['RSD','NO_RSD']:
        DRQ_filename = args.out_dir + '/master_picca_{}.fits'.format(RSD_option)
        catalog.write_DRQ(DRQ_filename,RSD_option,master_data,args.nside,overwrite=args.overwrite)

print('\nCreating the output file structure...')
#Make the new file structure
pixel_list = list(sorted(set(master_data['PIXNUM'])))
utils.make_file_structure(args.out_dir,pixel_list)
print(' -> Done!')
