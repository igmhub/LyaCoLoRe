#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from multiprocessing import Pool
import multiprocessing
import time
import argparse

from pyacolore import utils, catalog

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

parser.add_argument('--desi-footprint', action="store_true", default = False, required=False,
                    help = 'use the desi footprint (requires desimodel)')

parser.add_argument('--desi-footprint-plus', action="store_true", default = False, required=False,
                    help = 'use the desi footprint plus neighbouring pixels (requires desimodel)')

parser.add_argument('--downsampling', type = float, default = 1.0, required=False,
                    help = 'fraction by which to subsample the CoLoRe output')

parser.add_argument('--overwrite', action="store_true", default = False, required=False,
                    help = 'overwrite existing files')

################################################################################

args = parser.parse_args()

#Define global variables.

original_file_location = args.in_dir
new_base_file_location = args.out_dir
N_side = args.nside
min_catalog_z = args.min_cat_z
N_processes = args.nproc
parameter_filename = args.param_file
N_skewers = args.nskewers
add_picca_drqs = args.add_picca_drqs
desi_footprint = args.desi_footprint
desi_footprint_plus = args.desi_footprint_plus
if desi_footprint or desi_footprint_plus:
    desi_footprint = args.desi_footprint
    desi_footprint_plus = args.desi_footprint_plus
    try:
        from desimodel.footprint import tiles2pix
    except ModuleNotFoundError:
        raise InputError('Unable to use DESI footprint: desimodel is not installed.')
downsampling = args.downsampling
overwrite = args.overwrite

# TODO: print to confirm the arguments. e.g. "DLAs will be added"

#Check the value of N_side required is a power of 2.
if np.log2(N_side)-int(np.log2(N_side)) != 0:
    print('nside must be a power of 2!')
else:
    N_pix = 12*N_side**2

#Define the original file structure
original_filename_structure = 'out_srcs_s1_{}.fits' #file_number
file_numbers = list(range(0,32))
input_format = 'gaussian_colore'

#Set file structure
new_file_structure = '{}/{}/'               #pixel number//100, pixel number
new_filename_structure = '{}-{}-{}.fits'    #file type, nside, pixel number

#Get the simulation parameters from the parameter file.
simulation_parameters = utils.get_simulation_parameters(original_file_location,parameter_filename)

# TODO: Modify this to accomodate other density types.
#If density type is not lognormal, then crash.
if simulation_parameters['dens_type'] != 0:
    error('Density is not lognormal. Non-lognormal densities are not currently supported.')

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

#Define an error-tracking function.
def log_error(retval):
    print('Error:',retval)

################################################################################

"""
Produce the data required to make a master file using multiprocessing (one task per original file).
Join the outputs from the many processes together.
Save the master file, and a similarly structured file containing QSOs with 'bad coordinates'.
"""

print('\nWorking on master data...')
start = time.time()

#Choose the pixels we want.
if desi_footprint:
    pixels = tiles2pix(N_side)
elif desi_footprint_plus:
    pixels = tiles2pix(N_side)
    pixels = utils.add_pixel_neighbours(pixels)
else:
    pixels = None

#Define the process to make the master data.
def make_master_data(original_file_location,original_filename_structure,file_number,input_format,N_side,minimum_z=min_catalog_z):

    file_number, ID_data, cosmology, file_pixel_map_element, MOCKID_lookup_element = catalog.get_ID_data(original_file_location,original_filename_structure,file_number,input_format,N_side,minimum_z=min_catalog_z,pixels=pixels,downsampling=downsampling)

    return [file_number, ID_data, cosmology, file_pixel_map_element, MOCKID_lookup_element]

#Set up the multiprocessing pool parameters and make a list of tasks.
tasks = [(original_file_location,original_filename_structure,file_number,input_format,N_side) for file_number in file_numbers]

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = N_processes)
    results = []
    start_time = time.time()

    for task in tasks:
        pool.apply_async(make_master_data,task,callback=log_result,error_callback=log_error)

    pool.close()
    pool.join()

print('\nSaving the master files...')

#Join the multiprocessing results into 'master' and 'bad_coordinates' arrays.
master_data, bad_coordinates_data, cosmology_data, file_pixel_map, MOCKID_lookup = catalog.join_ID_data(results,N_side)

#Write master and bad coordinates files.
master_filename = new_base_file_location + '/master.fits'
catalog.write_ID(master_filename,master_data,cosmology_data,N_side,overwrite=overwrite)
print('\nMaster file contains {} objects.'.format(master_data.shape[0]))

if bad_coordinates_data.shape[0] > 0:
    bad_coordinates_filename = new_base_file_location + '/bad_coordinates.fits'
    catalog.write_ID(bad_coordinates_filename,bad_coordinates_data,cosmology_data,N_side,overwrite=overwrite)
    print('"bad coordinates" file contains {} objects.'.format(bad_coordinates_data.shape[0]))

#If desired, write the DRQ files for picca xcf to deal with.
if add_picca_drqs:
    print('\nMaster file contains {} objects.'.format(master_data.shape[0]))
    for RSD_option in ['RSD','NO_RSD']:
        DRQ_filename = new_base_file_location + '/master_picca_{}.fits'.format(RSD_option)
        catalog.write_DRQ(DRQ_filename,RSD_option,master_data,N_side,overwrite=overwrite)

print('Time to make master files: {:4.0f}s.'.format(time.time()-start))

#Make the new file structure
pixel_list = list(sorted(set(master_data['PIXNUM'])))
utils.make_file_structure(new_base_file_location,pixel_list)
