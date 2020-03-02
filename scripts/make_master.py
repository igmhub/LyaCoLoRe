#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from multiprocessing import Pool
import multiprocessing
import time
import argparse
import glob

from lyacolore import utils, catalog

################################################################################

#Script to produce a master file from CoLoRe's output files.

# TODO: Make 'input parameters' and 'output parameters' objects? Currently passing loads of arguments to multiprocessing functions which is messy
# TODO: option to reduce the number of skewers for test purposes?
# TODO: update the master file's cosmology once ssgf have been added
# TODO: use args to pass things like file structures around neatly?
# TODO: Get rid of the need to specify file numbers?

################################################################################

"""
Set up the file locations and filename structures.
Also define option preferences.
"""

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--in-dir', type = str, default = None, required=True,
                    help = 'input data directory')

parser.add_argument('--out-dir', type = str, default = None, required=True,
                    help = 'output data directory')

parser.add_argument('--file-format', type = str, default = 'colore', required=False,
                    choices=['colore'],
                    help = 'input file type')

parser.add_argument('--skewer-type', type = str, default = 'gaussian', required=False,
                    choices=['gaussian','density'],
                    help = 'type of skewer in input file')

parser.add_argument('--nproc', type = int, default = 1, required=False,
                    help = 'number of processes to use')

parser.add_argument('--nside', type = int, default = 16, required=False,
                    help = 'HEALPix nside for output files (must be 2^n)')

parser.add_argument('--min-cat-z', type = float, default = 1.8, required=False,
                    help = 'minimum z of objects in catalog')

parser.add_argument('--param-file', type = str, default = 'out_params.cfg', required=False,
                    help = 'output parameter file name')

parser.add_argument('--nskewers', type = int, default = None, required=False,
                    help = 'number of skewers to process')

parser.add_argument('--add-picca-drqs', action="store_true", default = False, required=False,
                    help = 'save picca format drq files')

parser.add_argument('--pixels', type = int, default = None, required=False,
                    help = 'which pixel numbers to work on', nargs='*')

parser.add_argument('--footprint', type = str, default = 'full_sky', required = False,
                    choices=['full_sky','desi','desi_pixel','desi_pixel_plus'],
                    help = 'name of footprint to use')

parser.add_argument('--downsampling', type = float, default = 1.0, required=False,
                    help = 'fraction by which to subsample the CoLoRe output')

parser.add_argument('--overwrite', action="store_true", default = False, required=False,
                    help = 'overwrite existing files')

parser.add_argument('--seed', type = int, default = 123, required=False,
                    help = 'specify seed to generate random numbers')

################################################################################

args = parser.parse_args()

# TODO: print to confirm the arguments. e.g. "DLAs will be added"

#Check the value of N_side required is a power of 2.
if np.log2(args.nside)-int(np.log2(args.nside)) != 0:
    print('nside must be a power of 2!')
else:
    N_pix = 12*args.nside**2

#Define the original file structure
input_filename_structure = 'out_srcs_s1_{}.fits' #file_number
input_files = glob.glob(args.in_dir+input_filename_structure.format('*'))
file_numbers = utils.get_file_numbers(args.in_dir,input_filename_structure,input_files)

#Set file structure
new_file_structure = '{}/{}/'               #pixel number//100, pixel number
new_filename_structure = '{}-{}-{}.fits'    #file type, nside, pixel number

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
def make_master_data(file_name,file_number,file_format,skewer_type,N_side,minimum_z=args.min_cat_z):

    seed = args.seed*10**5 + file_number
    file_number, ID_data, cosmology, file_pixel_map_element, MOCKID_lookup_element = catalog.get_ID_data(file_name,file_number,file_format,skewer_type,N_side,minimum_z=minimum_z,downsampling=args.downsampling,QSO_filter=QSO_filter,pixel_list=args.pixels,seed=seed)

    return [file_number, ID_data, cosmology, file_pixel_map_element, MOCKID_lookup_element]

#Set up the multiprocessing pool parameters and make a list of tasks.
tasks = [(input_files[i],file_number,args.file_format,args.skewer_type,args.nside) for i,file_number in enumerate(file_numbers)]

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
