#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import process_functions as functions
from multiprocessing import Pool
import multiprocessing
import sys
import time
import os

################################################################################

#Top level script to manage the conversion of CoLoRe output files into files more useful for analysis.
#These are produced on a per-HEALPix pixel basis, and currently include:
#Gaussian delta, density delta and flux delta files in picca format
#Gaussian and physical density files in CoLoRe format
#Transmission files

# TODO: Make 'input parameters' and 'output parameters' objects? Currently passing loads of arguments to multiprocessing functions which is messy

################################################################################

"""
Set up the file locations and filename structures.
Also define option preferences.
"""

#Define global variables.
lya = 1215.67

#Define the desired power of 2 for Nside for the output. This should be larger than that of the input.
N_side_pow2 = 3
N_side = 2**N_side_pow2
N_pix = 12*N_side**2

#Define the original file structure
original_file_location = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/output_G_hZ_4096_32_sr2.0_bm1/'
#original_file_location = '/Users/jfarr/Projects/test_data/output_G_hZ_4096_32_sr2.0_bm1/'
#original_file_location = '/Users/James/Projects/test_data/output_G_hZ_4096_32_sr2.0_bm1/'
original_filename_structure = 'out_srcs_s1_{}.fits' #file_number
file_numbers = list(range(0,1))
input_format = 'gaussian_colore'

#Set file structure
new_base_file_location = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/process_output_G_hZ_4096_32_sr2.0_bm1_nside{}_TEST/'.format(N_side)
#new_base_file_location = '/Users/jfarr/Projects/test_data/process_output_G_hZ_4096_32_sr2.0_bm1_nside8/'
#new_base_file_location = '/Users/James/Projects/test_data/test_adding_ssp/'
new_file_structure = '{}/{}/'               #pixel number//100, pixel number
new_filename_structure = '{}-{}-{}.fits'    #file type, nside, pixel number

#Choose options
minimum_catalog_z = 1.8
lambda_min = 3550 #A
zero_mean_delta = False
IVAR_cutoff = 1150 #A
final_cell_size = 0.25 #Mpc/h

# TODO: currently this only works for TAU= (A1*((1+z)/A2)**A3) * (density**alpha). Want it to also work for any A(z). Either import from a txt file or be able to point it to any function?
#Determine the different sets of parameters to test the density-flux conversion.
flux_parameters = [(0.374,4.0,5.1,1.0)] #A1,A2,A3,alpha

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

print('\nWorking on master data...')
start = time.time()

#Define the process to make the master data.
def make_master_data(original_file_location,original_filename_structure,file_number,input_format,N_side,minimum_z=minimum_catalog_z):

    file_number, ID_data, cosmology, file_pixel_map_element, MOCKID_lookup_element = functions.get_ID_data(original_file_location,original_filename_structure,file_number,input_format,N_side,minimum_z=minimum_catalog_z)

    return [file_number, ID_data, cosmology, file_pixel_map_element, MOCKID_lookup_element]

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

print('\nSaving the master files...')

#Join the multiprocessing results into 'master' and 'bad_coordinates' arrays.
master_data, bad_coordinates_data, cosmology_data, file_pixel_map, MOCKID_lookup = functions.join_ID_data(results,N_side)

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

#Make the new file structure
functions.make_file_structure(new_base_file_location,pixel_list)

print('\nWorking on per-HEALPix pixel initial Gaussian skewer files...')
start_time = time.time()

#Define the pixelisation process.
def pixelise_gaussian_skewers(pixel,original_file_location,original_filename_structure,input_format,shared_MOCKID_lookup,z_min,new_base_file_location,new_file_structure,N_side):
    #Define the save location for the pixel, according to the new file structure.
    location = new_base_file_location + new_file_structure.format(pixel//100,pixel)

    #Make file into an object
    pixel_object = functions.make_gaussian_pixel_object(pixel,original_file_location,original_filename_structure,input_format,shared_MOCKID_lookup,IVAR_cutoff=IVAR_cutoff)

    # TODO: These could be made beforehand and passed to the function? Or is there already enough being passed?
    #Make some useful headers
    header = fits.Header()
    header['NSIDE'] = N_side
    header['PIXNUM'] = pixel
    header['LYA'] = lya
    header['NQSO'] = pixel_object.N_qso

    # TODO: MISLEADING
    header['SIGMA_G'] = pixel_object.SIGMA_G

    #Gaussian CoLoRe
    filename = new_filename_structure.format('gaussian-colore',N_side,pixel)
    pixel_object.save_as_gaussian_colore(location,filename,header)

    #Don't need this here?
    """
    #Picca Gaussian
    filename = new_filename_structure.format('picca-gaussian',N_side,pixel)
    pixel_object.save_as_picca_gaussian(location,filename,header,zero_mean_delta=zero_mean_delta,lambda_min=lambda_min)
    """

    #Calculate the means of the pixel's gaussian skewers.
    #WARNING: this currently just uses all of the cells but this may be too slow once we've added small scale power?
    N, mean_DG, mean_DGS = functions.return_means(pixel_object.GAUSSIAN_DELTA_rows,pixel_object.IVAR_rows)
    means_data = [N,mean_DG,mean_DGS]

    return means_data

#Set up the multiprocessing pool parameters and make a list of tasks.
N_processes = int(sys.argv[1])
manager = multiprocessing.Manager()
shared_MOCKID_lookup = manager.dict(MOCKID_lookup)
tasks = [(pixel,original_file_location,original_filename_structure,input_format,shared_MOCKID_lookup,z_min,new_base_file_location,new_file_structure,N_side) for pixel in pixel_list]

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = N_processes)
    results = []
    start_time = time.time()

    for task in tasks:
        pool.apply_async(pixelise_gaussian_skewers,task,callback=log_result,error_callback=log_error)

    pool.close()
    pool.join()

print('\nTime to make Gaussian pixel files: {:4.0f}s.'.format(time.time()-start_time))
print(' ')

################################################################################

"""
To correctly calculate the physical fields, we must measure sigma from the Gaussian skewers.
"""

means_data_array = np.array(results)
N_total = np.sum(means_data_array[:,0])
gaussian_mean = (np.sum(means_data_array[:,1]*means_data_array[:,0]))/N_total
gaussian_variance = (np.sum(means_data_array[:,2]*means_data_array[:,0]))/N_total - gaussian_mean**2
measured_SIGMA_G = np.sqrt(gaussian_variance)

print('\nGaussian skewers have mean {:2.2f}, variance {:2.2f}.'.format(gaussian_mean,measured_SIGMA_G))

################################################################################

"""
We would like to add small scale flucatuations to the Gaussian field.
We must work out how much variance to add.
This is done by
 - computing the analytical P1D flux variance from Palanque-Delabrouille et al. (2013),
 - computing the Gaussian variance required to achieve this flux variance using Font-Ribera et al. (2012), and thus the extra variance our Gaussian skewers require
 - stretching the current skewers to achieve smaller cell sizes (NGP)
 - adding a random number to each cell with the appropriate statistics
"""
"""
#Calculate P1D variance
k = np.logspace(-5,5,10**6)
pk_z2 = aP1D.P1D_z_kms_PD2013(2.0,k)
W = np.sinc((k*0.25)/(2*np.pi))
sigma2_P1D = np.trapz((W**2)*pk_z2,k)
"""
#Work out sigma_G desired to achive the P1D sigma_F
# TODO: implement this from tune_flux.ipynb
desired_sigma_G = 4.0

#Determine the desired sigma_G by sampling
extra_sigma_G = np.sqrt(desired_sigma_G**2 - measured_SIGMA_G**2)

################################################################################

"""
We may now calculate the density and flux fields, and save the relevant files.
"""

print('\nWorking on per-HEALPix pixel final skewer files...')
start_time = time.time()

def produce_final_skewers(new_base_file_location,new_file_structure,new_filename_structure,pixel,N_side,input_format,zero_mean_delta,lambda_min,SIGMA_G):

    location = new_base_file_location + new_file_structure.format(pixel//100,pixel)

    #We work from the gaussian colore files made in 'pixelise gaussian skewers'
    gaussian_filename = new_filename_structure.format('gaussian-colore',N_side,pixel)

    #Make a pixel object from it
    pixel_object = functions.simulation_data.get_gaussian_skewers(location+gaussian_filename,None,input_format,SIGMA_G=measured_SIGMA_G,IVAR_cutoff=IVAR_cutoff)

    # TODO: These could be made beforehand and passed to the function? Or is there already enough being passed?
    #Make some useful headers
    header = fits.Header()
    header['NSIDE'] = N_side
    header['PIXNUM'] = pixel
    header['LYA'] = lya
    header['NQSO'] = pixel_object.N_qso

    #Add the physical density and flux skewers to the object.
    pixel_object.compute_physical_skewers() #Opportunity to change density type here
    pixel_object.compute_flux_skewers() #Opportunity to vary A and alpha here

    #lognorm CoLoRe
    filename = new_filename_structure.format('physical-colore',N_side,pixel)
    pixel_object.save_as_physical_colore(location,filename,header)

    #transmission
    filename = new_filename_structure.format('transmission',N_side,pixel)
    pixel_object.save_as_transmission(location,filename,header,lambda_min=lambda_min)

    #Add small scale power to the gaussian skewers:
    pixel_object.add_small_scale_gaussian_fluctuations(final_cell_size,extra_sigma_G,white_noise=True,lambda_min=0)

    #Picca Gaussian
    filename = new_filename_structure.format('picca-gaussian',N_side,pixel)
    pixel_object.save_as_picca_gaussian(location,filename,header,zero_mean_delta=zero_mean_delta,lambda_min=lambda_min,overwrite=True)

    #Picca density
    filename = new_filename_structure.format('picca-density',N_side,pixel)
    pixel_object.save_as_picca_density(location,filename,header,zero_mean_delta=zero_mean_delta,lambda_min=lambda_min)

    #picca flux
    filename = new_filename_structure.format('picca-flux',N_side,pixel)
    pixel_object.save_as_picca_flux(location,filename,header,lambda_min=lambda_min)

    return pixel

#define the tasks
tasks = [(new_base_file_location,new_file_structure,new_filename_structure,pixel,N_side,input_format,zero_mean_delta,lambda_min,measured_SIGMA_G) for pixel in pixel_list]

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = N_processes)
    results = []
    start_time = time.time()

    for task in tasks:
        pool.apply_async(produce_final_skewers,task,callback=log_result,error_callback=log_error)

    pool.close()
    pool.join()

print('\nTime to make physical pixel files: {:4.0f}s.'.format(time.time()-start_time))
print(' ')

################################################################################

# TODO: Update this
"""
Group the statistics calculated to get means and variances.
Save these into a fits file.
"""
"""
print('\nMaking statistics file...')
start_time = time.time()

#Use combine_means and means_to_statistics to calculate the mean and variance of the different quantities over all skewers.
means = functions.combine_means(results)
statistics = functions.means_to_statistics(means)

#Save the statistics data as a new fits file.
functions.write_statistics(new_base_file_location,N_side,statistics,cosmology_data)

print('\nTime to make statistics file: {:4.0f}s.'.format(time.time()-start_time))
print(' ')
"""

################################################################################

"""
Celebrate!
"""
