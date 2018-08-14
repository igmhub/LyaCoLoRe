#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from multiprocessing import Pool
import multiprocessing
import sys
import time
import os
import argparse

from pyacolore import general, independent, stats, convert, pixelise, input, DLA, RSD

################################################################################

#Script to make a transmission file for a given pixel, given a

# TODO: Make 'input parameters' and 'output parameters' objects? Currently passing loads of arguments to multiprocessing functions which is messy
# TODO: option to reduce the number of skewers for test purposes?
# TODO: update the master file's cosmology once ssgf have been added
# TODO: update handling of sigma_G as it's a bit misleading atm
# TODO: use args to pass things like file structures around neatly?
# TODO: Get rid of the need to specify file numbers?
# TODO: Issue with if we don't have consecutive pixel numbers?

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

parser.add_argument('--master-dir', type = str, default = None, required=False,
                    help = 'directory containing the master file')

parser.add_argument('--nproc', type = int, default = 1, required=False,
                    help = 'number of processes to use')

parser.add_argument('--nside', type = int, default = 16, required=False,
                    help = 'HEALPix nside for output files (must be 2^n)')

parser.add_argument('--pixels', type = int, default = None, required=False,
                    help = 'which pixel numbers to work on', nargs='*')

parser.add_argument('--IVAR-cut', type = float, default = 1150., required=False,
                    help = 'maximum rest frame lambda for IVAR=1 (Å)')

parser.add_argument('--cell-size', type = float, default = 0.25, required=False,
                    help = 'size in Mpc/h of output cells')

parser.add_argument('--lambda-min', type = float, default = 3550., required=False,
                    help = 'minimum lambda in picca skewers (Å)')

parser.add_argument('--min-cat-z', type = float, default = 1.8, required=False,
                    help = 'minimum z of objects in catalog')

parser.add_argument('--param-file', type = str, default = 'out_params.cfg', required=False,
                    help = 'output parameter file name')

parser.add_argument('--tuning-file', type = str, default = 'input_files/tune_small_scale_fluctuations.fits', required=False,
                    help = 'file name for data about tuning sigma_G/alpha')

parser.add_argument('--add-DLAs', action="store_true", default = True, required=False,
                    help = 'add DLAs to the transmission file')

parser.add_argument('--add-RSDs', action="store_true", default = False, required=False,
                    help = 'add linear RSDs to the transmission file')

parser.add_argument('--add-Lyb', action="store_true", default = False, required=False,
                    help = 'add Lyman-beta absorption to TRANSMISSION HDU')

parser.add_argument('--add-metals', action="store_true", default = False, required=False,
                    help = 'add METALS HDU with metal absorption')

parser.add_argument('--include-thermal-effects', action="store_true", default = False, required=False,
                    help = 'add thermal RSDs to the transmission file')

parser.add_argument('--retune-small-scale-fluctuations', action="store_true", default = False, required=False,
                    help = 'recalculate the values of sigma_G and alpha needed')

parser.add_argument('--transmission-only', action="store_true", default = False, required=False,
                    help = 'save only the transmission file')

parser.add_argument('--nskewers', type = int, default = None, required=False,
                    help = 'number of skewers to process')

parser.add_argument('--seed', type = int, default = 123, required=False,
                    help = 'specify seed to generate random numbers')

################################################################################

args = parser.parse_args()

#Define global variables.
lya = 1215.67

original_file_location = args.in_dir
new_base_file_location = args.out_dir
master_location = args.master_dir
if not master_location:
    master_location = new_base_file_location
N_side = args.nside
pixels = args.pixels
if not pixels:
    pixels = list(range(12*N_side**2))
min_catalog_z = args.min_cat_z
lambda_min = args.lambda_min
zero_mean_delta = False
IVAR_cutoff = args.IVAR_cut
final_cell_size = args.cell_size
N_processes = args.nproc
parameter_filename = args.param_file
add_DLAs = args.add_DLAs
add_RSDs = args.add_RSDs
add_Lyb = args.add_Lyb
add_metals = args.add_metals
include_thermal_effects = args.include_thermal_effects
retune_small_scale_fluctuations = args.retune_small_scale_fluctuations
tuning_file = args.tuning_file
transmission_only = args.transmission_only
N_skewers = args.nskewers
global_seed = args.seed

# TODO: print to confirm the arguments. e.g. "DLAs will be added"

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

#Calculate the minimum value of z that we are interested in.
#i.e. the z value for which lambda_min cooresponds to the lya wavelength.
z_min = lambda_min/lya - 1

#Get the simulation parameters from the parameter file.
simulation_parameters = general.get_simulation_parameters(original_file_location,parameter_filename)

# TODO: Modify this to accomodate other density types.
#If density type is not lognormal, then crash.
if simulation_parameters['dens_type'] != 0:
    error('Density is not lognormal. Non-lognormal densities are not currently supported.')

################################################################################

"""
Construct the MOCKID_lookup from the master file.
"""
# TODO: potential issue with differnt values of nside being used in make_master.py
master = fits.open(master_location+'/master.fits')
master_data = master[1].data
master.close()

#Make a MOCKID lookup.
pixel_list = list(sorted(set([pixel for pixel in master_data['PIXNUM'] if pixel in pixels])))
MOCKID_lookup = {}
for pixel in pixel_list:
    #pixel_indices = [i for i in range(len(master_data['PIXNUM'])) if master_data['PIXNUM'][i]==pixel]
    pixel_indices = (master_data['PIXNUM']==pixel)
    pixel_MOCKIDs = master_data['MOCKID'][pixel_indices]
    pixel_file_number_list = list(sorted(set(master_data['FILENUM'][pixel_indices])))
    for file_number in pixel_file_number_list:
        pixel_file_indices = ((master_data['FILENUM'][pixel_indices])==file_number)
        #MOCKID_list = [master_data['MOCKID'][i] for i in range(len(master_data['PIXNUM'])) if master_data['PIXNUM'][i]==pixel and master_data['FILENUM'][i]==file_number]
        pixel_file_MOCKIDs = pixel_MOCKIDs[pixel_file_indices]
        MOCKID_lookup = {**MOCKID_lookup,**{(file_number,pixel):list(pixel_file_MOCKIDs)}}

################################################################################

"""
Define the multiprocessing tracking functions
"""

#Define a progress-tracking function.
def log_result(retval):

    results.append(retval)
    N_complete = len(results)
    N_tasks = len(tasks)

    general.progress_bar(N_complete,N_tasks,start_time)

#Define an error-tracking function.
def log_error(retval):
    print('Error:',retval)

################################################################################

print('\nWorking on per-HEALPix pixel initial Gaussian skewer files...')
start_time = time.time()

#Define the pixelisation process.
def pixelise_gaussian_skewers(pixel,original_file_location,original_filename_structure,input_format,MOCKID_lookup,z_min,new_base_file_location,new_file_structure,N_side):
    #Define the save location for the pixel, according to the new file structure.
    location = new_base_file_location + '/' + new_file_structure.format(pixel//100,pixel)

    #Make file into an object
    pixel_object = pixelise.make_gaussian_pixel_object(pixel,original_file_location,original_filename_structure,input_format,MOCKID_lookup,IVAR_cutoff=IVAR_cutoff)

    # TODO: These could be made beforehand and passed to the function? Or is there already enough being passed?
    #Make some useful headers
    header = fits.Header()
    header['HPXNSIDE'] = N_side
    header['HPXPIXEL'] = pixel
    header['HPXNEST'] = True
    header['LYA'] = lya

    # TODO: MISLEADING
    header['SIGMA_G'] = pixel_object.SIGMA_G

    #Gaussian CoLoRe
    filename = new_filename_structure.format('gaussian-colore',N_side,pixel)
    pixel_object.save_as_gaussian_colore(location,filename,header)

    #Calculate the means of the pixel's gaussian skewers.
    #WARNING: this currently just uses all of the cells but this may be too slow once we've added small scale power?
    N, mean_DG, mean_DGS = stats.return_means(pixel_object.GAUSSIAN_DELTA_rows,pixel_object.IVAR_rows)
    means_data = [N,mean_DG,mean_DGS]

    return means_data

#Set up the multiprocessing pool parameters and make a list of tasks.
#what's the sharing doing here?
#manager = multiprocessing.Manager()
#shared_MOCKID_lookup = manager.dict(MOCKID_lookup)
tasks = [(pixel,original_file_location,original_filename_structure,input_format,MOCKID_lookup,z_min,new_base_file_location,new_file_structure,N_side) for pixel in pixel_list]

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = N_processes)
    results = []
    start_time = time.time()

    for task in tasks:
        pool.apply_async(pixelise_gaussian_skewers,task,callback=log_result,error_callback=log_error)

    pool.close()
    pool.join()

print('\nTime to make Gaussian pixel files: {:4.0f}s.\n'.format(time.time()-start_time))

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
If desired, we can recalculate the tuning parameters.
Otherwise, we just load values from file.
"""

#Work out sigma_G desired to achive the P1D sigma_dF
beta = 1.65

if retune_small_scale_fluctuations == True:

    # TODO: this needs to be written.
    import tune_flux_parameters
    tune_small_scale_fluctuations = tune_flux_parameters.tune()

else:
    #Otherwise, load the data from the fits file that has been pre-computed.
    print('\nLoading how much extra power to add from file...')
    h = fits.open(tuning_file)
    tune_small_scale_fluctuations = h['DATA'].data
    h.close()
    print('Process complete!')

tuning_z_values = tune_small_scale_fluctuations['z']
desired_sigma_G_values = tune_small_scale_fluctuations['sigma_G']
desired_mean_F = tune_small_scale_fluctuations['mean_F']
tuning_alphas = tune_small_scale_fluctuations['alpha']

#Determine the desired sigma_G by sampling
extra_sigma_G_values = np.sqrt(desired_sigma_G_values**2 - measured_SIGMA_G**2)

################################################################################

"""
We may now calculate the density and flux fields, and save the relevant files.
"""

print('\nWorking on per-HEALPix pixel final skewer files...')
start_time = time.time()

def produce_final_skewers(new_base_file_location,new_file_structure,new_filename_structure,pixel,N_side,input_format,zero_mean_delta,lambda_min,measured_SIGMA_G):
    location = new_base_file_location + '/' + new_file_structure.format(pixel//100,pixel)
    mean_F_data = np.array(list(zip(tuning_z_values,desired_mean_F)))

    #We work from the gaussian colore files made in 'pixelise gaussian skewers'.
    gaussian_filename = new_filename_structure.format('gaussian-colore',N_side,pixel)

    #Make a pixel object from it.
    pixel_object = pixelise.SimulationData.get_gaussian_skewers_object(location+gaussian_filename,None,input_format,SIGMA_G=measured_SIGMA_G,IVAR_cutoff=IVAR_cutoff)

    if add_Lyb:
        pixel_object.setup_Lyb_absorber()

    if add_metals:
        pixel_object.setup_metal_absorbers()

    #Make some useful headers
    header = fits.Header()
    header['HPXNSIDE'] = N_side
    header['HPXPIXEL'] = pixel
    header['HPXNEST'] = True
    header['LYA'] = lya
    header['SIGMA_G'] = measured_SIGMA_G

    if transmission_only == False:
        #lognorm CoLoRe
        pixel_object.compute_physical_skewers()
        filename = new_filename_structure.format('physical-colore',N_side,pixel)
        pixel_object.save_as_physical_colore(location,filename,header)

    #Trim the skewers (remove low lambda cells)
    pixel_object.trim_skewers(lambda_min,min_catalog_z,extra_cells=1)

    #Exit now if no skewers are left.
    if pixel_object.N_qso == 0:
        print('\nwarning: no objects left in pixel {} after trimming.'.format(pixel))
        return pixel

    #Get seed to generate random numbers for this particular pixel
    #seed = 10**(len(str(12*N_side**2))) + pixel + global_seed
    seed = int(str(N_side) + str(pixel)) + global_seed

    #Add small scale power to the gaussian skewers:
    generator = np.random.RandomState(seed)
    new_cosmology = pixel_object.add_small_scale_gaussian_fluctuations(final_cell_size,tuning_z_values,extra_sigma_G_values,generator,white_noise=False,lambda_min=lambda_min,IVAR_cutoff=IVAR_cutoff)
    #new_cosmology = []

    #Remove the 'SIGMA_G' header as SIGMA_G now varies with z, so can't be stored in a header.
    del header['SIGMA_G']
    pixel_object.SIGMA_G = np.interp(pixel_object.Z,tuning_z_values,desired_sigma_G_values)

    #Add a table with DLAs in to the pixel object.
    # TODO: in future, we want DLAs all the way down to z=0.
    #That means we need to store skewers all the way down to z=0.
    #Not possible atm as we'd run out of memory, but can be done once running on >1 node.
    if add_DLAs:
        pixel_object.add_DLA_table(seed)

    #Add physical skewers to the object.
    pixel_object.compute_physical_skewers()

    #Add tau skewers to the object, starting with Lyman-alpha
    alphas=np.interp(pixel_object.Z,tuning_z_values,tuning_alphas)
    pixel_object.compute_all_tau_skewers(alphas,beta)

    if transmission_only == False:
        #Picca Gaussian
        filename = new_filename_structure.format('picca-gaussian',N_side,pixel)
        pixel_object.save_as_picca_gaussian(location,filename,header)

        #Picca density
        filename = new_filename_structure.format('picca-density',N_side,pixel)
        pixel_object.save_as_picca_density(location,filename,header)

        #picca flux
        filename = new_filename_structure.format('picca-flux-noRSD',N_side,pixel)
        pixel_object.save_as_picca_flux(location,filename,header,mean_F_data=mean_F_data)

    #Add thermal RSDs to the tau skewers.
    #Add RSDs from the velocity skewers provided by CoLoRe.
    if add_RSDs == True:
        pixel_object.add_all_RSDs(alphas,beta,thermal=include_thermal_effects)

    #transmission
    filename = new_filename_structure.format('transmission',N_side,pixel)
    pixel_object.save_as_transmission(location,filename,header)

    if transmission_only == False:
        #Picca Gaussian
        #filename = new_filename_structure.format('picca-gaussian',N_side,pixel)
        #pixel_object.save_as_picca_gaussian(location,filename,header)

        #Picca density
        #filename = new_filename_structure.format('picca-density',N_side,pixel)
        #pixel_object.save_as_picca_density(location,filename,header)

        #picca flux
        filename = new_filename_structure.format('picca-flux',N_side,pixel)
        pixel_object.save_as_picca_flux(location,filename,header,mean_F_data=mean_F_data)
    else:
        #If transmission_only is not False, remove the gaussian-colore file.
        os.remove(location+gaussian_filename)

    means = pixel_object.get_means()

    return [new_cosmology,means]

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

print('\nTime to make physical pixel files: {:4.0f}s.\n'.format(time.time()-start_time))

################################################################################
"""
Having added small scale power, we must add a new HDU to the master file's cosmology.
"""

print('Updating master file\'s cosmology...')
#First check that the new cosmologies are all the same.
# TODO: some kind of system to check consistency here?
new_cosmology = results[0][0]
master_filename = master_location + '/master.fits'

#Reorganise the data.
master = fits.open(master_filename)
try:

    test = master[3].data
    master.close()

except IndexError:

    master_catalog = master[1].data
    master_colore_cosmology = master[2].data
    master_new_cosmology = new_cosmology

    #Make an appropriate header.
    header = fits.Header()
    header['NSIDE'] = N_side

    #Make the data into tables.
    hdu_ID = fits.BinTableHDU.from_columns(master_catalog,header=header,name='CATALOG')
    hdu_cosmology_colore = fits.BinTableHDU.from_columns(master_colore_cosmology,header=header,name='COSMO_COL')
    hdu_cosmology_expanded = fits.BinTableHDU.from_columns(master_new_cosmology,header=header,name='COSMO_EXP')

    #Make a primary HDU.
    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)

    #Make the .fits file.
    hdulist = fits.HDUList([prihdu,hdu_ID,hdu_cosmology_colore,hdu_cosmology_expanded])
    hdulist.writeto(master_filename,overwrite=True)
    hdulist.close()

print('Process complete!\n')

################################################################################

"""
Group the statistics calculated to get means and variances.
Save these into a fits file.
"""
"""
print('\nMaking statistics file...')
start_time = time.time()

#Use combine_means and means_to_statistics to calculate the mean and variance of the different quantities over all skewers.
means_list = []
for result in results:
    means_list += [result[1]]
means = stats.combine_means(means_list)
statistics = stats.means_to_statistics(means)

#Save the statistics data as a new fits file.
functions.write_statistics(new_base_file_location,N_side,statistics,new_cosmology)

print('\nTime to make statistics file: {:4.0f}s.\n'.format(time.time()-start_time))
"""
################################################################################

"""
Celebrate!
"""
