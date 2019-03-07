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

from pyacolore import utils, independent, stats, convert, simulation_data, DLA, RSD, tuning

################################################################################

#Script to make a transmission file for a given pixel, given a

# TODO: Make 'input parameters' and 'output parameters' objects? Currently passing loads of arguments to multiprocessing functions which is messy
# TODO: option to reduce the number of skewers for test purposes?
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

parser.add_argument('--tuning-file', type = str, default = 'input_files/tuning_data_151118.fits', required=False,
                    help = 'file name for data about tuning sigma_G/alpha')

parser.add_argument('--add-small-scale-fluctuations', action="store_true", default = False, required=False,
                    help = 'add small scale fluctuations to the Gaussian skewers')

parser.add_argument('--add-DLAs', action="store_true", default = False, required=False,
                    help = 'add DLAs to the transmission file')

parser.add_argument('--DLA-bias', type = float, default = 2., required=False,
                    help = 'bias of DLAs')

parser.add_argument('--DLA-bias-method', type = str, default = 'b_const', required=False,
                    help = 'either \"b_const\" or \"bD_const\"')

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

parser.add_argument('--overwrite', action="store_true", default = False, required=False,
                    help = 'overwrite existing files')

parser.add_argument('--add-QSO-RSDs', action="store_true", default = False, required=False,
                    help = 'add QSO RSDs to the transmission file')

parser.add_argument('--smoothing-R-kms', type = float, default = 25.0, required=False,
                    help = 'size in km/s of extra power smoothing radius')

# TODO: this is now defunct.
parser.add_argument('--fit-function-to-tuning-data', action="store_true", default = False, required=False,
                    help = 'fit a function of the form A0 * (z^A1) + A2 to the tuning data')

################################################################################

print('setup arguments from parser')

args = parser.parse_args()

#Define global variables.
lya = utils.lya_rest

base_in_dir = args.in_dir
base_out_dir = args.out_dir
if args.master_dir:
    master_file = args.master_dir+'/master.fits'
else:
    master_file = base_out_dir+'/master.fits'
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
add_ssf = args.add_small_scale_fluctuations
add_DLAs = args.add_DLAs
dla_bias = args.DLA_bias
dla_bias_method = args.DLA_bias_method
add_RSDs = args.add_RSDs
add_Lyb = args.add_Lyb
add_metals = args.add_metals
include_thermal_effects = args.include_thermal_effects
retune_small_scale_fluctuations = args.retune_small_scale_fluctuations
tuning_file = args.tuning_file
transmission_only = args.transmission_only
N_skewers = args.nskewers
global_seed = args.seed
fit_function_to_tuning_data = args.fit_function_to_tuning_data
overwrite = args.overwrite
add_QSO_RSDs = args.add_QSO_RSDs
R_kms = args.smoothing_R_kms

# TODO: print to confirm the arguments. e.g. "DLAs will be added"

if np.log2(N_side)-int(np.log2(N_side)) != 0:
    print('nside must be a power of 2!')
else:
    N_pix = 12*N_side**2

#colore skewers filename (except number that will be added later)
colore_base_filename = base_in_dir+'/out_srcs_s1_'

#Calculate the minimum value of z that we are interested in.
#i.e. the z value for which lambda_min cooresponds to the lya wavelength.
z_min = lambda_min/lya - 1

#Get the simulation parameters from the parameter file.
simulation_parameters = utils.get_simulation_parameters(base_in_dir,parameter_filename)

# TODO: Modify this to accomodate other density types.
#If density type is not lognormal, then crash.
if simulation_parameters['dens_type'] != 0:
    error('Density is not lognormal. Non-lognormal densities are not currently supported.')
input_format='gaussian_colore'

################################################################################

"""
Construct the MOCKID_lookup from the master file.
"""
# TODO: potential issue with differnt values of nside being used in make_master.py
master = fits.open(master_file)
master_data = master[1].data
master.close()

#Make a MOCKID lookup.
master_data_pixel_set = set(master_data['PIXNUM'])
pixels_set = set(pixels)
pixel_list = list(sorted(master_data_pixel_set.intersection(pixels_set)))

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

    utils.progress_bar(N_complete,N_tasks,start_time)

#Define an error-tracking function.
def log_error(retval):
    print('Error:',retval)

################################################################################

print('\nWorking on per-HEALPix pixel initial Gaussian skewer files...')
start_time = time.time()

#Define the pixelisation process.
def pixelise_gaussian_skewers(pixel,colore_base_filename,z_min,base_out_dir,N_side):

    #Define the output directory the pixel, according to the new file structure.
    location = utils.get_dir_name(base_out_dir,pixel)

    #at some point we might want to read physical density
    input_format='gaussian_colore'

    #Make file into an object
    pixel_object = simulation_data.make_gaussian_pixel_object(pixel,colore_base_filename,input_format,shared_MOCKID_lookup,IVAR_cutoff=IVAR_cutoff)

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
    filename = utils.get_file_name(location,'gaussian-colore',N_side,pixel)
    pixel_object.save_as_colore('gaussian',filename,header,overwrite=overwrite)

    #Calculate the means of the pixel's gaussian skewers.
    N = np.sum(pixel_object.IVAR_rows,axis=0)
    mean_DG = np.average(pixel_object.get_mean_quantity('gaussian',power=1),weights=N)
    mean_DGS = np.average(pixel_object.get_mean_quantity('gaussian',power=2),weights=N)
    N = np.sum(N)

    return (N,mean_DG,mean_DGS)

#Set up the multiprocessing pool parameters and make a list of tasks.
#what's the sharing doing here?
manager = multiprocessing.Manager()
shared_MOCKID_lookup = manager.dict(MOCKID_lookup)
tasks = [(pixel,colore_base_filename,z_min,base_out_dir,N_side) for pixel in pixel_list]

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

gaussian_mean = np.average(means_data_array[:,1],weights=means_data_array[:,0])
gaussian_variance = np.average(means_data_array[:,2],weights=means_data_array[:,0]) - gaussian_mean**2
measured_SIGMA_G = np.sqrt(gaussian_variance)

print('\nGaussian skewers have mean {:2.2f}, variance {:2.2f}.'.format(gaussian_mean,measured_SIGMA_G))
print('\nModifying header showing sigma_G in Gaussian CoLoRe files...')

def modify_header(pixel):
    location = utils.get_dir_name(base_out_dir,pixel)
    filename = utils.get_file_name(location,'gaussian-colore',N_side,pixel)
    h = fits.open(filename)
    for HDU in h[1:]:
        HDU.header['SIGMA_G'] = measured_SIGMA_G
    h.writeto(filename,overwrite=True)
    h.close()
    return

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = N_processes)
    results = []
    start_time = time.time()

    for pixel in pixel_list:
        pool.apply_async(modify_header,(pixel,),callback=log_result,error_callback=log_error)

    pool.close()
    pool.join()

################################################################################

"""
We would like to add small scale flucatuations to the Gaussian field.
We load values of the parameters from file.
"""

h = fits.open(tuning_file)
n = h[1].header['n']
k1 = h[1].header['k1']
tuning_z_values = h[1].data['z']
tuning_alphas = h[1].data['alpha']
tuning_betas = h[1].data['beta']
tuning_sigma_Gs = h[1].data['sigma_G']
h.close()

################################################################################

"""
We may now do the main work of LyaCoLoRe. This includes:
 - add extra small scale power
 - convert from the gaussian field to the lognormal field, then tau
 - add metals
 - add RSDs
 - add DLAs
 - convert from tau to flux
 - save the transmission files
We also save picca format delta files for running correlation function tests.
Deltas are caclulated using the mean quantity in each pixel.
They are renormalised using the global mean in 'make_summaries'
"""

print('\nWorking on per-HEALPix pixel final skewer files...')
start_time = time.time()

def produce_final_skewers(base_out_dir,pixel,N_side,zero_mean_delta,lambda_min,measured_SIGMA_G,n,k1):

    t = time.time()

    #Define a random seed for use in this pixel.
    seed = int(pixel * 10**5 + global_seed)

    #We work from the gaussian colore files made in 'pixelise gaussian skewers'.
    location = utils.get_dir_name(base_out_dir,pixel)
    gaussian_filename = utils.get_file_name(location,'gaussian-colore',N_side,pixel)

    #Make a pixel object from it.
    file_number = None
    pixel_object = simulation_data.SimulationData.get_gaussian_skewers_object(gaussian_filename,file_number,input_format,SIGMA_G=measured_SIGMA_G,IVAR_cutoff=IVAR_cutoff)

    #print('{:3.2f} checkpoint object'.format(time.time()-t)); t = time.time()

    #Add Lyb and metal absorbers if needed.
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

    #Save CoLoRe format files.
    if transmission_only == False:
        #lognorm CoLoRe
        pixel_object.compute_physical_skewers()
        filename = utils.get_file_name(location,'physical-colore',N_side,pixel)
        pixel_object.save_as_colore('density',filename,header,overwrite=overwrite)

    #Trim the skewers (remove low lambda cells). Exit if no QSOs are left.
    #We don't cut too tightly on the low lambda to allow for RSDs.
    lambda_buffer = 100. #A
    pixel_object.trim_skewers(lambda_min-lambda_buffer,min_catalog_z,extra_cells=1)
    if pixel_object.N_qso == 0:
        print('\nwarning: no objects left in pixel {} after trimming.'.format(pixel))
        return pixel

    #Save picca format files without adding small scale power.
    if transmission_only == False:
        filename = utils.get_file_name(location,'picca-gaussian-colorecell',N_side,pixel)
        pixel_object.save_as_picca_delta('gaussian',filename,header,overwrite=overwrite)

        filename = utils.get_file_name(location,'picca-density-colorecell',N_side,pixel)
        pixel_object.save_as_picca_delta('density',filename,header,overwrite=overwrite)

    #print('{:3.2f} checkpoint colore files'.format(time.time()-t)); t = time.time()

    #Add a table with DLAs in to the pixel object.
    # TODO: in future, we want DLAs all the way down to z=0.
    #That means we need to store skewers all the way down to z=0.
    #May need to adjust how many nodes are used when running.
    if add_DLAs:
        pixel_object.add_DLA_table(seed,dla_bias=dla_bias,method=dla_bias_method)

    #print('{:3.2f} checkpoint DLAs'.format(time.time()-t)); t = time.time()

    #Add small scale power to the gaussian skewers:
    if add_ssf:
        generator = np.random.RandomState(seed)
        pixel_object.add_small_scale_gaussian_fluctuations(final_cell_size,tuning_z_values,tuning_sigma_Gs,generator,white_noise=False,lambda_min=lambda_min,IVAR_cutoff=IVAR_cutoff,n=n,k1=k1,R_kms=R_kms)

        #Remove the 'SIGMA_G' header as SIGMA_G now varies with z, so can't be stored in a header.
        del header['SIGMA_G']
        sigma_G = np.sqrt(tuning_sigma_Gs**2 + measured_SIGMA_G**2)
        pixel_object.SIGMA_G = np.exp(np.interp(np.log(pixel_object.Z),np.log(tuning_z_values),np.log(sigma_G)))

    #print('{:3.2f} checkpoint SSF'.format(time.time()-t)); t = time.time()

    #Recompute physical skewers.
    pixel_object.compute_physical_skewers()

    #Add tau skewers to the object, starting with Lyman-alpha
    alphas = np.exp(np.interp(np.log(pixel_object.Z),np.log(tuning_z_values),np.log(tuning_alphas)))
    betas = np.exp(np.interp(np.log(pixel_object.Z),np.log(tuning_z_values),np.log(tuning_betas)))
    sigma_Gs = np.exp(np.interp(np.log(pixel_object.Z),np.log(tuning_z_values),np.log(tuning_sigma_Gs)))
    pixel_object.compute_all_tau_skewers(alphas,betas)

    if transmission_only == False:

        #Picca Gaussian, small cells
        filename = utils.get_file_name(location,'picca-gaussian',N_side,pixel)
        pixel_object.save_as_picca_delta('gaussian',filename,header,overwrite=overwrite,add_QSO_RSDs=add_QSO_RSDs)

        #Picca density
        filename = utils.get_file_name(location,'picca-density',N_side,pixel)
        pixel_object.save_as_picca_delta('density',filename,header,overwrite=overwrite,add_QSO_RSDs=add_QSO_RSDs)

        #Picca tau
        filename = utils.get_file_name(location,'picca-tau-noRSD-notnorm',N_side,pixel)
        pixel_object.save_as_picca_delta('tau',filename,header,notnorm=True,overwrite=overwrite,add_QSO_RSDs=False)

        #Picca flux
        filename = utils.get_file_name(location,'picca-flux-noRSD-notnorm',N_side,pixel)
        pixel_object.save_as_picca_delta('flux',filename,header,notnorm=True,overwrite=overwrite,add_QSO_RSDs=False)

    #Save the no RSD statistics file for this pixel.
    filename = 'statistics-noRSD-16-{}.fits'.format(pixel)
    statistics = pixel_object.save_statistics(location,filename,overwrite=overwrite)

    #print('{:3.2f} checkpoint noRSD files'.format(time.time()-t)); t = time.time()

    #Add RSDs from the velocity skewers provided by CoLoRe.
    if add_RSDs == True:
        pixel_object.add_all_RSDs(thermal=include_thermal_effects)

    #print('{:3.2f} checkpoint RSDs'.format(time.time()-t)); t = time.time()

    #Trim the skewers (remove low lambda cells). Exit if no QSOs are left.
    #We now cut hard at lambda min as RSDs have been implemented.
    pixel_object.trim_skewers(lambda_min,min_catalog_z,extra_cells=1)
    if pixel_object.N_qso == 0:
        print('\nwarning: no objects left in pixel {} after trimming.'.format(pixel))
        return pixel

    #Make a variable containing the new cosmology data.
    new_cosmology = pixel_object.return_cosmology()

    #transmission
    filename = utils.get_file_name(location,'transmission',N_side,pixel)
    pixel_object.save_as_transmission(filename,header,overwrite=overwrite)

    if transmission_only == False:
        #Picca tau
        filename = utils.get_file_name(location,'picca-tau-notnorm',N_side,pixel)
        pixel_object.save_as_picca_delta('tau',filename,header,notnorm=True,overwrite=overwrite,add_QSO_RSDs=add_QSO_RSDs)

        #Picca flux
        filename = utils.get_file_name(location,'picca-flux-notnorm',N_side,pixel)
        pixel_object.save_as_picca_delta('flux',filename,header,notnorm=True,overwrite=overwrite,add_QSO_RSDs=add_QSO_RSDs)
    else:
        #If transmission_only is not False, remove the gaussian-colore file.
        os.remove(gaussian_filename)

    #Save the final statistics file for this pixel.
    filename = 'statistics-16-{}.fits'.format(pixel)
    statistics = pixel_object.save_statistics(location,filename,overwrite=overwrite)

    #print('{:3.2f} checkpoint RSD files'.format(time.time()-t)); t = time.time()

    return new_cosmology

#define the tasks
tasks = [(base_out_dir,pixel,N_side,zero_mean_delta,lambda_min,measured_SIGMA_G,n,k1) for pixel in pixel_list]

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
PROBABLY COULD MOVE THIS TO make_summaries
Having added small scale power, we must add a new HDU to the master file's cosmology.
"""

print('Updating master file\'s cosmology...')
#First check that the new cosmologies are all the same.
# TODO: some kind of system to check consistency here?
new_cosmology = results[0]

#Reorganise the data.
master = fits.open(master_file)
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
    hdulist.writeto(master_file,overwrite=True)
    hdulist.close()

print('Process complete!\n')

################################################################################

"""
Celebrate!
"""
