#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import process_functions as functions
import matplotlib.pyplot as plt
from multiprocessing import Pool
import multiprocessing
import sys
import time
import os
import argparse

################################################################################

#Top level script to manage the conversion of CoLoRe output files into files more useful for analysis.
#These are produced on a per-HEALPix pixel basis, and currently include:
#Gaussian delta, density delta and flux delta files in picca format
#Gaussian and physical density files in CoLoRe format
#Transmission files

# TODO: Make 'input parameters' and 'output parameters' objects? Currently passing loads of arguments to multiprocessing functions which is messy
# TODO: option to reduce the number of skewers for test purposes?
# TODO: update the master file's cosmology once ssgf have been added
# TODO: update handling of sigma_G as it's a bit misleading atm
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

parser.add_argument('--add-RSDs', action="store_true", default = True, required=False,
                    help = 'add RSDs to the transmission file')

parser.add_argument('--retune-small-scale-fluctuations', action="store_true", default = False, required=False,
                    help = 'recalculate the values of sigma_G and alpha needed')

parser.add_argument('--transmission-only', action="store_true", default = False, required=False,
                    help = 'save only the transmission file')

parser.add_argument('--nskewers', type = int, default = None, required=False,
                    help = 'number of skewers to process')

################################################################################

args = parser.parse_args()

#Define global variables.
lya = 1215.67

original_file_location = args.in_dir
new_base_file_location = args.out_dir
N_side = args.nside
min_catalog_z = args.min_cat_z
lambda_min = args.lambda_min
zero_mean_delta = False
IVAR_cutoff = args.IVAR_cut
final_cell_size = args.cell_size
N_processes = args.nproc
parameter_filename = args.param_file
add_DLAs = args.add_DLAs
add_RSDs = args.add_RSDs
retune_small_scale_fluctuations = args.retune_small_scale_fluctuations
tuning_file = args.tuning_file
transmission_only = args.transmission_only
N_skewers = args.nskewers

# TODO: print to confirm the arguments. e.g. "DLAs will be added"

if np.log2(N_side)-int(np.log2(N_side)) != 0:
    print('nside must be a power of 2!')
else:
    N_pix = 12*N_side**2

#Define the original file structure
original_filename_structure = 'N1000_out_srcs_s1_{}.fits' #file_number
file_numbers = list(range(0,1))
input_format = 'gaussian_colore'

#Set file structure
new_file_structure = '{}/{}/'               #pixel number//100, pixel number
new_filename_structure = '{}-{}-{}.fits'    #file type, nside, pixel number

#Calculate the minimum value of z that we are interested in.
#i.e. the z value for which lambda_min cooresponds to the lya wavelength.
z_min = lambda_min/lya - 1

#Get the simulation parameters from the parameter file.
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
    print(' -> current progress: {} {:4d} of {:4d} complete ({:3.0%}), {:4.0f}s elapsed, ~{:5.0f}s remaining'.format(progress_bar,N_complete,N_tasks,N_complete/N_tasks,time_elapsed,estimated_time_remaining),end="\r")

    if len(results) == len(tasks):
        print('\nProcess complete!')

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
def make_master_data(original_file_location,original_filename_structure,file_number,input_format,N_side,minimum_z=min_catalog_z):

    file_number, ID_data, cosmology, file_pixel_map_element, MOCKID_lookup_element = functions.get_ID_data(original_file_location,original_filename_structure,file_number,input_format,N_side,minimum_z=min_catalog_z)

    return [file_number, ID_data, cosmology, file_pixel_map_element, MOCKID_lookup_element]

#Set up the multiprocessing pool parameters and make a list of tasks.
#N_processes = int(sys.argv[1])
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
master_filename = new_base_file_location + '/master.fits'
functions.write_ID(master_filename,master_data,cosmology_data,N_side)
print('\nMaster file contains {} objects.'.format(master_data.shape[0]))

if bad_coordinates_data.shape[0] > 0:
    bad_coordinates_filename = new_base_file_location + '/bad_coordinates.fits'
    functions.write_ID(bad_coordinates_filename,bad_coordinates_data,cosmology_data,N_side)
    print('"bad coordinates" file contains {} objects.'.format(bad_coordinates_data.shape[0]))

#Write the DRQ files for picca xcf to deal with.
for RSD_option in ['RSD','NO_RSD']:
    DRQ_filename = new_base_file_location + '/master_picca_{}.fits'.format(RSD_option)
    functions.write_DRQ(DRQ_filename,RSD_option,master_data,N_side)

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
    location = new_base_file_location + '/' + new_file_structure.format(pixel//100,pixel)

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

    #Calculate the means of the pixel's gaussian skewers.
    #WARNING: this currently just uses all of the cells but this may be too slow once we've added small scale power?
    N, mean_DG, mean_DGS = functions.return_means(pixel_object.GAUSSIAN_DELTA_rows,pixel_object.IVAR_rows)
    means_data = [N,mean_DG,mean_DGS]

    return means_data

#Set up the multiprocessing pool parameters and make a list of tasks.
#N_processes = int(sys.argv[1])
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
We must work out how much variance to add.
This is done by
 - computing the analytical P1D flux variance from Palanque-Delabrouille et al. (2013),
 - computing the Gaussian variance required to achieve this flux variance using Font-Ribera et al. (2012), and thus the extra variance our Gaussian skewers require
 - stretching the current skewers to achieve smaller cell sizes (NGP)
 - adding a random number to each cell with the appropriate statistics
"""

#Work out sigma_G desired to achive the P1D sigma_dF
tuning_z_values = np.linspace(0,4.0,128)
beta = 1.65

if retune_small_scale_fluctuations == True:

    print('\nCalculating how much extra power to add...')

    D_values = np.interp(tuning_z_values,cosmology_data['Z'],cosmology_data['D'])
    sigma_G_tolerance = 0.0001

    def tune_sigma_G(z,D,l_hMpc,beta,Om):

        sigma_dF_needed = functions.get_sigma_dF_P1D(z,l_hMpc=l_hMpc,Om=Om)
        mean_F_needed = functions.get_mean_F_model(z)

        alpha,sigma_G,mean_F,sigma_dF = functions.find_sigma_G(mean_F_needed,sigma_dF_needed,beta,D,tolerance=sigma_G_tolerance)

        return (z,alpha,sigma_G,mean_F,sigma_dF,mean_F_needed,sigma_dF_needed)

    tasks = [(z,np.interp(z,cosmology_data['Z'],cosmology_data['D']),final_cell_size,beta,simulation_parameters['omega_M']) for z in tuning_z_values]

    if __name__ == '__main__':
        pool = Pool(processes = N_processes)
        results = []
        start_time = time.time()

        for task in tasks:
            pool.apply_async(tune_sigma_G,task,callback=log_result,error_callback=log_error)

        pool.close()
        pool.join()

    dtype = [('z', 'f4'), ('alpha', 'f4'), ('sigma_G', 'f4'), ('mean_F', 'f4'), ('sigma_dF', 'f4'), ('mean_F_needed', 'f4'), ('sigma_dF_needed', 'f4')]
    tune_small_scale_fluctuations = np.array(results,dtype=dtype)
    tune_small_scale_fluctuations = np.sort(tune_small_scale_fluctuations,order=['z'])

    plt.plot(tune_small_scale_fluctuations['z'],tune_small_scale_fluctuations['alpha'],label='alpha')
    plt.plot(tune_small_scale_fluctuations['z'],tune_small_scale_fluctuations['sigma_G'],label='sigma_G')
    plt.plot(tune_small_scale_fluctuations['z'],tune_small_scale_fluctuations['mean_F'],label='mean_F')
    plt.plot(tune_small_scale_fluctuations['z'],tune_small_scale_fluctuations['sigma_dF'],label='sigma_dF')
    plt.grid()
    plt.legend()
    plt.savefig('tune_flux_values_tol{}_n{}.pdf'.format(sigma_G_tolerance,tuning_z_values.shape[0]))
    plt.show()
    """
    plt.plot(tune_small_scale_fluctuations['z'],tune_small_scale_fluctuations['mean_F']/tune_small_scale_fluctuations['mean_F_needed'] - 1,label='mean_F error')
    plt.plot(tune_small_scale_fluctuations['z'],tune_small_scale_fluctuations['sigma_dF']/tune_small_scale_fluctuations['sigma_dF_needed'] - 1,label='sigma_dF error')
    plt.grid()
    plt.legend()
    plt.savefig('tune_flux_values_tol{}_n{}_Ferrors.pdf'.format(sigma_G_tolerance,tuning_z_values.shape[0]))
    #plt.show()
    """
    header = fits.Header()
    header['beta'] = beta
    header['l_hMpc'] = final_cell_size

    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)
    cols_DATA = fits.ColDefs(tune_small_scale_fluctuations)
    hdu_DATA = fits.BinTableHDU.from_columns(cols_DATA,header=header,name='DATA')

    hdulist = fits.HDUList([prihdu, hdu_DATA])
    hdulist.writeto('input_files/tune_small_scale_fluctuations_n{}.fits'.format(tuning_z_values.shape[0]))
    hdulist.close()

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
alphas = tune_small_scale_fluctuations['alpha']

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
    pixel_object = functions.simulation_data.get_gaussian_skewers_object(location+gaussian_filename,None,input_format,SIGMA_G=measured_SIGMA_G,IVAR_cutoff=IVAR_cutoff)

    #Make some useful headers
    header = fits.Header()
    header['NSIDE'] = N_side
    header['PIXNUM'] = pixel
    header['LYA'] = lya
    header['NQSO'] = pixel_object.N_qso
    header['NESTED'] = True
    header['SIGMA_G'] = measured_SIGMA_G

    #Add the physical density and flux skewers to the object.
    pixel_object.compute_physical_skewers()
    pixel_object.compute_flux_skewers(np.interp(pixel_object.Z,tuning_z_values,alphas),beta)

    if transmission_only == False:
        #lognorm CoLoRe
        filename = new_filename_structure.format('physical-colore',N_side,pixel)
        pixel_object.save_as_physical_colore(location,filename,header)

    #Trim the skewers (remove low lambda cells)
    pixel_object.trim_skewers(lambda_min,min_catalog_z,extra_cells=1)

    #Exit now if no skewers are left.
    if pixel_object.N_qso == 0:
        print('\nwarning: no objects left in pixel {} after trimming.'.format(pixel))
        return pixel

    #Save picca files for testing RSDs
    filename = new_filename_structure.format('picca-gaussian-noRSD',N_side,pixel)
    pixel_object.save_as_picca_gaussian(location,filename,header)
    filename = new_filename_structure.format('picca-flux-noRSD',N_side,pixel)
    pixel_object.save_as_picca_flux(location,filename,header,mean_F_data=mean_F_data)

    #Add RSDs from the velocity skewers provided by CoLoRe.
    if add_RSDs == True:
        pixel_object.add_linear_RSDs()

    #If there are already physical/flux skewers, recalculate them.
    if pixel_object.DENSITY_DELTA_rows is not None:
        pixel_object.compute_physical_skewers()
    if pixel_object.F_rows is not None:
        pixel_object.compute_flux_skewers(np.interp(pixel_object.Z,tuning_z_values,alphas),beta)

    #Save picca files for testing RSDs
    filename = new_filename_structure.format('picca-gaussian-RSD',N_side,pixel)
    pixel_object.save_as_picca_gaussian(location,filename,header)
    filename = new_filename_structure.format('picca-flux-RSD',N_side,pixel)
    pixel_object.save_as_picca_flux(location,filename,header,mean_F_data=mean_F_data)

    #Add small scale power to the gaussian skewers:
    new_cosmology = pixel_object.add_small_scale_gaussian_fluctuations(final_cell_size,tuning_z_values,extra_sigma_G_values,white_noise=True,lambda_min=lambda_min,IVAR_cutoff=IVAR_cutoff)

    #Remove the 'SIGMA_G' header as SIGMA_G now varies with z.
    del header['SIGMA_G']

    #If there are already physical/flux skewers, recalculate them.
    if pixel_object.DENSITY_DELTA_rows is not None:
        pixel_object.compute_physical_skewers()
    if pixel_object.F_rows is not None:
        pixel_object.compute_flux_skewers(np.interp(pixel_object.Z,tuning_z_values,alphas),beta)

    #Add a table with DLAs in to the pixel object.
    # TODO: in future, we want DLAs all the way down to z=0.
    #That means we need to store skewers all the way down to z=0.
    #Not possible atm as we'd run out of memory, but can be done once running on >1 node.
    if add_DLAs:
        pixel_object.add_DLA_table()

    #transmission
    filename = new_filename_structure.format('transmission',N_side,pixel)
    pixel_object.save_as_transmission(location,filename,header)

    if transmission_only == False:
        #Picca Gaussian
        filename = new_filename_structure.format('picca-gaussian',N_side,pixel)
        pixel_object.save_as_picca_gaussian(location,filename,header)

        #Picca density
        filename = new_filename_structure.format('picca-density',N_side,pixel)
        pixel_object.save_as_picca_density(location,filename,header)

        #picca flux
        filename = new_filename_structure.format('picca-flux',N_side,pixel)
        pixel_object.save_as_picca_flux(location,filename,header,mean_F_data=mean_F_data)

        #picca velocity
        filename = new_filename_structure.format('picca-velocity',N_side,pixel)
        pixel_object.save_as_picca_velocity(location,filename,header)
    else:
        #If transmission_only is not False, remove the gaussian-colore file
        os.remove(location+gaussian_filename)

    means = pixel_object.get_means(lambda_min=lambda_min)

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

#Reorganise the data.
master = fits.open(master_filename)
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

print('\nMaking statistics file...')
start_time = time.time()

#Use combine_means and means_to_statistics to calculate the mean and variance of the different quantities over all skewers.
means_list = list(np.array(results)[:,1])
means = functions.combine_means(means_list)
statistics = functions.means_to_statistics(means)

#Save the statistics data as a new fits file.
functions.write_statistics(new_base_file_location,N_side,statistics,master_new_cosmology)

print('\nTime to make statistics file: {:4.0f}s\n.'.format(time.time()-start_time))


################################################################################

"""
We may want to renormalise some of the deltas.
Here we do so.
"""
"""
file_to_renormalise_names = ['picca-gaussian-RSD']
file_to_renormalise_types = ['picca-gaussian']

tasks = []
for i,filename in enumerate(file_to_renormalise_names):
    file_type = file_to_renormalise_types[i]
    tasks += [(filename,file_type,pixel) for pixel in pixel_list]

def renormalise_delta_file(filename,file_type,pixel):
    location = new_base_file_location + '/' + new_file_structure.format(pixel//100,pixel)
    filename
    pixel_object = get_gaussian_skewers_object(location+filename,None,input_format,SIGMA_G=measured_SIGMA_G,IVAR_cutoff=IVAR_cutoff)
"""




################################################################################


"""
Celebrate!
"""
