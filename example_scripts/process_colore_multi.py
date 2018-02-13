#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import process_functions as functions
from multiprocessing import Pool
import sys
import time
import os

################################################################################

"""
Set up the file locations and filename structures.
Also define option preferences.
"""

#Top level script to manage the conversion of CoLoRe output files into files more useful for analysis.
#These are produced on a per-HEALPix pixel basis.

lya = 1215.67

#Define the desired power of 2 for Nside for the output. This should be larger than that of the input.
N_side_pow2 = 3
N_side = 2**N_side_pow2
N_pix = 12*N_side**2

#Define the original file structure
original_file_location = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_revamp/output_4096_32/'
original_file_location = '/Users/jfarr/Projects/repixelise/test_input/'
original_filename_structure = 'N100000_out_srcs_s1_{}.fits' #file_number
file_numbers = list(range(15,16))
input_format = 'physical_colore'

#Set file structure
new_base_file_location = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_revamp/test/lya1100/'
new_base_file_location = '/Users/jfarr/Projects/repixelise/test_output/test_multi/'
new_file_structure = '{}/{}/'               #pixel number//100, pixel number
new_filename_structure = '{}-{}-{}.fits'    #file type, nside, pixel number

#Choose options
lambda_min = 3550
zero_mean_delta = False
IVAR_cutoff = 1150

#Calculate the minimum value of z that we are interested in.
#i.e. the z value for which lambda_min cooresponds to the lya wavelength.
z_min = lambda_min/lya - 1

#Get the simulation parameters from the parameter file.
parameter_filename = 'out_params.cfg'
simulation_parameters = functions.get_simulation_parameters(original_file_location,parameter_filename)

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
    N_tasks_digits = int(np.log10(N_tasks)) + 1

    N_chunks = 20
    N_chunks_complete = int((N_complete*N_chunks) // (N_tasks))
    block_char = '-'
    progress_bar = '|' + block_char*N_chunks_complete + ' '*(N_chunks-N_chunks_complete) + '|'

    current_time = time.time()
    time_elapsed = current_time - start_time
    estimated_time_remaining = (time_elapsed)*((N_tasks-N_complete)/N_complete)
    print(' -> current progress: {} {:4d} of {:4d} complete ({:3.0%}), {:4.0f}s elapsed, ~{:4.0f}s remaining'.format(progress_bar,N_complete,N_tasks,N_complete/N_tasks,time_elapsed,estimated_time_remaining),end="\r")

    if len(results) == len(tasks):
        print('\n\nProcess complete!')

#Define an error-tracking function.
def log_error(retval):
    print('Error:',retval)

################################################################################

"""
Produce the data required to make a master file using multiprocessing (one task per original file).
Join the outputs from the many processes together.
Save the master file, and a similarly structured file containing QSOs with 'bad coordinates'.
"""

print('\nWorking on master file...')
start = time.time()

#Define the process to make the master data.
def make_master_data(original_file_location,original_filename_structure,file_number,input_format,N_side):

    ID_data, cosmology = functions.get_ID_data(original_file_location,original_filename_structure,file_number,input_format,N_side)

    return [ID_data, cosmology]

#Set up the multiprocessing pool parameters and make a list of tasks.
N_processes = int(sys.argv[1])
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

#Join the multiprocessing results into 'master' and 'bad_coordinates' arrays.
master_data, bad_coordinates_data, cosmology_data = functions.join_ID_data(results)

#Make a list of the pixels that the files cover.
pixel_list = list(sorted(set(master_data['PIXNUM'])))
file_number_list = list(sorted(set(file_numbers)))

#Write master and bad coordinates files.
master_filename = new_base_file_location + '/' + 'nside_{}_'.format(N_side) + 'master.fits'
functions.write_ID(master_filename,master_data,cosmology_data,N_side)

bad_coordinates_filename = new_base_file_location + '/' + 'nside_{}_'.format(N_side) + 'bad_coordinates.fits'
functions.write_ID(bad_coordinates_filename,bad_coordinates_data,cosmology_data,N_side)

#Write the DRQ files for picca xcf to deal with.
for RSD_option in ['RSD','NO_RSD']:
    DRQ_filename = new_base_file_location + '/nside_{}_master_picca_{}.fits'.format(N_side,RSD_option)
    functions.write_DRQ(DRQ_filename,RSD_option,master_data,N_side)

print('\nMaster file contains {} objects.\n"bad coordinates" file contains {} objects.'.format(master_data.shape[0],bad_coordinates_data.shape[0]))
print('Time to make master file: {:4.0f}s.'.format(time.time()-start))

################################################################################

"""
Make the new file structure.
For each pixel, run 'pixelise' to make a pixel object, and save it in the required formats.
Done using a multiprocessing pool, with one task per pixel.
"""

print('\nWorking on per-HEALPix pixel skewer files...')
start = time.time()

#Make the new file structure
functions.make_file_structure(new_base_file_location,pixel_list)

#Define the pixelisation process.
def pixelise(pixel,original_file_location,original_filename_structure,input_format,master_data,pixel_list,file_number_list,z_min,new_base_file_location,new_file_structure,N_side):

    #Define the save location for the pixel, according to the new file structure.
    location = new_base_file_location + new_file_structure.format(pixel//100,pixel)

    #Make file into an object
    pixel_object = functions.make_pixel_object(pixel,original_file_location,original_filename_structure,input_format,master_data,pixel_list,file_number_list,IVAR_cutoff=IVAR_cutoff)

    #Make some useful headers
    header = fits.Header()
    header['NSIDE'] = N_side
    header['PIXNUM'] = pixel
    header['LYA'] = lya
    header['NQSO'] = pixel_object.N_qso
    header['SIGMA_G'] = pixel_object.SIGMA_G

    #Gaussian CoLoRe
    filename = new_filename_structure.format('gaussian-colore',N_side,pixel)
    pixel_object.save_as_gaussian_colore(location,filename,header)

    #Picca Gaussian
    filename = new_filename_structure.format('picca-gaussian',N_side,pixel)
    pixel_object.save_as_picca_gaussian(location,filename,header,zero_mean_delta=zero_mean_delta,lambda_min=lambda_min)

    #lognorm CoLoRe
    filename = new_filename_structure.format('physical-colore',N_side,pixel)
    pixel_object.save_as_physical_colore(location,filename,header)

    #Picca density
    filename = new_filename_structure.format('picca-density',N_side,pixel)
    pixel_object.save_as_picca_density(location,filename,header,zero_mean_delta=zero_mean_delta,lambda_min=lambda_min)

    #transmission
    filename = new_filename_structure.format('transmission',N_side,pixel)
    pixel_object.save_as_transmission(location,filename,header,lambda_min=lambda_min)

    #picca flux
    filename = new_filename_structure.format('picca-flux',N_side,pixel)
    pixel_object.save_as_picca_flux(location,filename,header,lambda_min=lambda_min)

    return pixel

#Set up the multiprocessing pool parameters and make a list of tasks.
N_processes = int(sys.argv[1])
tasks = [(pixel,original_file_location,original_filename_structure,input_format,master_data,pixel_list,file_number_list,z_min,new_base_file_location,new_file_structure,N_side) for pixel in pixel_list]

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = N_processes)
    results = []
    start_time = time.time()

    for task in tasks:
        pool.apply_async(pixelise,task,callback=log_result,error_callback=log_error)

    pool.close()
    pool.join()

################################################################################

"""
Celebrate!
"""

print('\nTime to make pixel files: {:4.0f}s.'.format(time.time()-start))
print(' ')
