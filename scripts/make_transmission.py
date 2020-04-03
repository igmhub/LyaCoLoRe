#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from multiprocessing import Pool
import multiprocessing
from scipy.interpolate import interp1d
import sys
import time
import os
import argparse

from lyacolore import utils, independent, stats, convert, simulation_data, DLA, RSD, tuning

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

parser.add_argument('--file-format', type = str, default = 'colore', required=False,
                    choices=['colore'],
                    help = 'input file type')

parser.add_argument('--skewer-type', type = str, default = 'gaussian', required=False,
                    choices=['gaussian','density'],
                    help = 'type of skewer in input file')

parser.add_argument('--master-dir', type = str, default = None, required=False,
                    help = 'directory containing the master file')

parser.add_argument('--nproc', type = int, default = 1, required=False,
                    help = 'number of processes to use')

parser.add_argument('--nside', type = int, default = 16, required=False,
                    help = 'HEALPix nside for output files (must be 2^n)')

parser.add_argument('--pixels', type = int, default = None, required=False,
                    help = 'which pixel numbers to work on', nargs='*')

parser.add_argument('--IVAR-cut', type = float, default = 1200., required=False,
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

parser.add_argument('--DLA-bias-evol', type = str, default = 'b_const', required=False,
                    choices=['b_const','bD_const'],
                    help = 'choose DLA bias evolution with redshift')

parser.add_argument('--DLA-bias-method', type = str, default = 'global', required=False,
                    choices=['global','sample'],
                    help = 'choose whether the DLA bias is determined by the global or sample value of sigma_G')

parser.add_argument('--add-RSDs', action="store_true", default = False, required=False,
                    help = 'add linear RSDs to the transmission file')

parser.add_argument('--add-Lyb', action="store_true", default = False, required=False,
                    help = 'add Lyman-beta absorption to TRANSMISSION HDU')

parser.add_argument('--add-metals', action="store_true", default = False, required=False,
                    help = 'include metal absorbers in the transmission files')

parser.add_argument('--picca-all-absorbers', action="store_true", default = False, required=False,
                    help = 'combine all absorbers in the picca tau and flux files')

parser.add_argument('--include-thermal-effects', action="store_true", default = False, required=False,
                    help = 'add thermal RSDs to the transmission file')

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

parser.add_argument('--transmission-lambda-min', type = float, default = 3550., required=False,
                    help = 'minimum wavelength stored in the transmission files')

parser.add_argument('--transmission-lambda-max', type = float, default = 6500., required=False,
                    help = 'maximum wavelength stored in the transmission files')

parser.add_argument('--transmission-delta-lambda', type = float, default = 0.2, required=False,
                    help = 'pixel size of transmission files wavelength grid')

parser.add_argument('--transmission-format', type = str, default = "final", required=False,
                    choices=['develop','final','single_HDU'],
                    help = 'format of transmission files')

parser.add_argument('--compress', action="store_true", default = False, required=False,
                    help = 'compress output files to .fits.gz')

################################################################################

print('setup arguments from parser')

args = parser.parse_args()

#Define global variables.
lya = utils.lya_rest

if args.master_dir:
    master_file = args.master_dir+'/master.fits'
else:
    master_file = args.out_dir+'/master.fits'
if not args.pixels:
    args.pixels = list(range(12*args.nside**2))

zero_mean_delta = False

# TODO: print to confirm the arguments. e.g. "DLAs will be added"

if np.log2(args.nside)-int(np.log2(args.nside)) != 0:
    print('nside must be a power of 2!')
else:
    N_pix = 12*args.nside**2

#colore skewers filename (except number that will be added later)
colore_base_filename = args.in_dir+'/out_srcs_s1_'

#Calculate the minimum value of z that we are interested in.
#i.e. the z value for which lambda_min cooresponds to the lya wavelength.
z_min = args.lambda_min/lya - 1
small = 10**-10

#Get the simulation parameters from the parameter file.
simulation_parameters = utils.get_simulation_parameters(args.in_dir,args.param_file)

#If we have density input skewers and want to add DLAs, then raise an error:
#this functionality is not yet implemented.
if (args.skewer_type=='density') & args.add_DLAs:
    raise ValueError('Adding DLAs from density input skewers is not possible yet!')

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
pixels_set = set(args.pixels)
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

print('\nWorking on per-HEALPix pixel initial skewer files...')
start_time = time.time()

#Define the pixelisation process.
def pixelise_colore_output(pixel,colore_base_filename,z_min,out_dir,N_side):

    #Define the output directory the pixel, according to the new file structure.
    location = utils.get_dir_name(out_dir,pixel)

    #Make file into an object
    pixel_object = simulation_data.make_pixel_object(pixel,colore_base_filename,args.file_format,args.skewer_type,shared_MOCKID_lookup,IVAR_cutoff=args.IVAR_cut)

    # TODO: These could be made beforehand and passed to the function? Or is there already enough being passed?
    #Make some useful headers
    header = fits.Header()
    header['HPXNSIDE'] = N_side
    header['HPXPIXEL'] = pixel
    header['HPXNEST'] = True
    header['LYA'] = lya

    ## Save the pixelised colore file.
    filename = utils.get_file_name(location,'{}-colore'.format(args.skewer_type),N_side,pixel)
    pixel_object.save_as_colore(args.skewer_type,filename,header,overwrite=args.overwrite,compress=args.compress)

    if args.skewer_type == 'gaussian':
        pixel_object.compute_SIGMA_G(type='single_value',lr_max=args.IVAR_cut)
        header['SIGMA_G'] = pixel_object.SIGMA_G
        N = np.sum(pixel_object.IVAR_rows.astype('int'))

        return (N,pixel_object.SIGMA_G)

    else:

        return

#Set up the multiprocessing pool parameters and make a list of tasks.
#what's the sharing doing here?
manager = multiprocessing.Manager()
shared_MOCKID_lookup = manager.dict(MOCKID_lookup)
tasks = [(pixel,colore_base_filename,z_min,args.out_dir,args.nside) for pixel in pixel_list]

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = args.nproc)
    results = []
    start_time = time.time()

    for task in tasks:
        pool.apply_async(pixelise_colore_output,task,callback=log_result,error_callback=log_error)

    pool.close()
    pool.join()

print('\nTime to make pixel files: {:4.0f}s.\n'.format(time.time()-start_time))

################################################################################

"""
To correctly calculate the physical fields when using Gaussian input skewers, we
must measure sigma from all skewers.
Here, we combine the results from each pixel, to compute an overall value.
"""

if args.skewer_type == 'gaussian':
    N_values = np.array([r[0] for r in results])
    sg_values = np.array([r[1] for r in results])
    SIGMA_G_global = np.sqrt(np.sum((sg_values**2)*N_values)/np.sum(N_values))

    print('\nGaussian skewers have mean sigma {:2.4f}.'.format(SIGMA_G_global))
    print('\nModifying header showing sigma_G in Gaussian CoLoRe files...')

    def modify_header(pixel):
        location = utils.get_dir_name(args.out_dir,pixel)
        filename = utils.get_file_name(location,'gaussian-colore',args.nside,pixel,compressed=args.compress)
        h = fits.open(filename)
        for HDU in h[1:]:
            HDU.header['SIGMA_G'] = SIGMA_G_global
        h.writeto(filename,overwrite=True)
        h.close()
        return

    #Run the multiprocessing pool
    if __name__ == '__main__':
        pool = Pool(processes = args.nproc)
        results = []
        start_time = time.time()

        for pixel in pixel_list:
            pool.apply_async(modify_header,(pixel,),callback=log_result,error_callback=log_error)

        pool.close()
        pool.join()


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
Deltas are normalised using the global mean in 'make_summaries'
"""

print('\nWorking on per-HEALPix pixel final skewer files...')
start_time = time.time()

def produce_final_skewers(base_out_dir,pixel,N_side,zero_mean_delta,lambda_min,tuning_file):

    t = time.time()

    # Define a random seed for use in this pixel.
    seed = int(pixel * 10**5 + args.seed)

    #We work from the gaussian colore files made in 'pixelise gaussian skewers'.
    location = utils.get_dir_name(base_out_dir,pixel)
    gaussian_filename = utils.get_file_name(location,'{}-colore'.format(args.skewer_type),N_side,pixel,compressed=args.compress)

    # Make a pixel object from it.
    file_number = None
    pixel_object = simulation_data.SimulationData.get_skewers_object(gaussian_filename,file_number,args.file_format,args.skewer_type,IVAR_cutoff=args.IVAR_cut)
    if args.skewer_type == 'gaussian':
        pixel_object.SIGMA_G = SIGMA_G_global

    # Make a transformation object and add it to the pixel object.
    pixel_object.add_transformation_from_file(tuning_file)

    #Scale the velocities.
    pixel_object.scale_velocities(use_transformation=True)

    #print('{:3.2f} checkpoint object'.format(time.time()-t)); t = time.time()

    #Add Lyb and metal absorbers if needed.
    if args.add_Lyb:
        pixel_object.setup_Lyb_absorber()
    if args.add_metals:
        pixel_object.setup_metal_absorbers()

    #Make some useful headers
    header = fits.Header()
    header['HPXNSIDE'] = N_side
    header['HPXPIXEL'] = pixel
    header['HPXNEST'] = True
    header['LYA'] = lya
    if args.skewer_type == 'gaussian':
        header['SIGMA_G'] = pixel_object.SIGMA_G

    #Save CoLoRe format files.
    if args.transmission_only == False:
        if args.skewer_type == 'gaussian':
            pixel_object.compute_physical_skewers()
        filename = utils.get_file_name(location,'density-colore',N_side,pixel)
        pixel_object.save_as_colore('density',filename,header,overwrite=args.overwrite,compress=args.compress)

    #Trim the skewers (remove low lambda cells). Exit if no QSOs are left.
    #We don't cut too tightly on the low lambda to allow for RSDs.
    lambda_buffer = 100. #A
    pixel_object.trim_skewers(lambda_min-lambda_buffer,args.min_cat_z,extra_cells=1)
    if pixel_object.N_qso == 0:
        print('\nwarning: no objects left in pixel {} after trimming.'.format(pixel))
        return pixel

    #Save picca format files without adding small scale power.
    if args.transmission_only == False:
        if args.skewer_type == 'gaussian':
            filename = utils.get_file_name(location,'picca-gaussian-colorecell',N_side,pixel)
            pixel_object.save_as_picca_delta('gaussian',filename,header,overwrite=args.overwrite,add_QSO_RSDs=args.add_QSO_RSDs,compress=args.compress)

        filename = utils.get_file_name(location,'picca-density-colorecell',N_side,pixel)
        pixel_object.save_as_picca_delta('density',filename,header,overwrite=args.overwrite,add_QSO_RSDs=args.add_QSO_RSDs,compress=args.compress)

    #print('{:3.2f} checkpoint colore files'.format(time.time()-t)); t = time.time()

    #Add a table with DLAs in to the pixel object.
    # TODO: in future, we want DLAs all the way down to z=0.
    #That means we need to store skewers all the way down to z=0.
    #May need to adjust how many nodes are used when running.

    if args.add_DLAs:
        pixel_object.add_DLA_table(seed,dla_bias=args.DLA_bias,evol=args.DLA_bias_evol,method=args.DLA_bias_method)

    #print('{:3.2f} checkpoint DLAs'.format(time.time()-t)); t = time.time()

    #Add small scale power to the gaussian skewers:
    if args.add_small_scale_fluctuations:
        generator = np.random.RandomState(seed)
        pixel_object.add_small_scale_fluctuations(args.cell_size,generator,white_noise=False,lambda_min=lambda_min,IVAR_cutoff=args.IVAR_cut,use_transformation=True)

        if args.skewer_type == 'gaussian':
            #Remove the 'SIGMA_G' header as SIGMA_G now varies with z, so can't be stored in a header.
            del header['SIGMA_G']

    #print('{:3.2f} checkpoint SSF'.format(time.time()-t)); t = time.time()

    #Recompute physical skewers, and then the tau skewers.
    if args.skewer_type == 'gaussian':
        pixel_object.compute_physical_skewers()
    pixel_object.compute_all_tau_skewers()

    if args.transmission_only == False:

        if args.skewer_type == 'gaussian':
            #Picca Gaussian, small cells
            filename = utils.get_file_name(location,'picca-gaussian',N_side,pixel)
            pixel_object.save_as_picca_delta('gaussian',filename,header,overwrite=args.overwrite,add_QSO_RSDs=args.add_QSO_RSDs,compress=args.compress)

        #Picca density
        filename = utils.get_file_name(location,'picca-density',N_side,pixel)
        pixel_object.save_as_picca_delta('density',filename,header,overwrite=args.overwrite,add_QSO_RSDs=args.add_QSO_RSDs,compress=args.compress)

        #Picca tau
        filename = utils.get_file_name(location,'picca-tau-noRSD-notnorm',N_side,pixel)
        pixel_object.save_as_picca_delta('tau',filename,header,notnorm=True,overwrite=args.overwrite,add_QSO_RSDs=args.add_QSO_RSDs,compress=args.compress,all_absorbers=args.picca_all_absorbers)

        #Picca flux
        filename = utils.get_file_name(location,'picca-flux-noRSD-notnorm',N_side,pixel)
        pixel_object.save_as_picca_delta('flux',filename,header,notnorm=True,overwrite=args.overwrite,add_QSO_RSDs=args.add_QSO_RSDs,compress=args.compress,all_absorbers=args.picca_all_absorbers)

        """
        ## Disable this for the moment.
        #Save the no RSD statistics file for this pixel.
        filename = utils.get_file_name(location,'statistics-noRSD',N_side,pixel)
        statistics = pixel_object.save_statistics(filename,overwrite=args.overwrite,compress=args.compress,all_absorbers=args.picca_all_absorbers)
        """

    #print('{:3.2f} checkpoint noRSD files'.format(time.time()-t)); t = time.time()

    #Add RSDs from the velocity skewers provided by CoLoRe.
    if args.add_RSDs == True:
        pixel_object.add_all_RSDs(thermal=args.include_thermal_effects)

    #print('{:3.2f} checkpoint RSDs'.format(time.time()-t)); t = time.time()

    #Trim the skewers (remove low lambda cells). Exit if no QSOs are left.
    #We now cut hard at lambda min as RSDs have been implemented.
    pixel_object.trim_skewers(lambda_min,args.min_cat_z,extra_cells=1)
    if pixel_object.N_qso == 0:
        print('\nwarning: no objects left in pixel {} after trimming.'.format(pixel))
        return pixel

    #Make a variable containing the new cosmology data.
    new_cosmology = pixel_object.return_cosmology()

    #Save the transmission file.
    filename = utils.get_file_name(location,'transmission',N_side,pixel)
    pixel_object.save_as_transmission(filename,header,overwrite=args.overwrite,wave_min=args.transmission_lambda_min,wave_max=args.transmission_lambda_max,wave_step=args.transmission_delta_lambda,fmt=args.transmission_format,add_QSO_RSDs=args.add_QSO_RSDs,compress=args.compress)

    if args.transmission_only == False and args.add_RSDs == True:
        #Picca tau
        filename = utils.get_file_name(location,'picca-tau-notnorm',N_side,pixel)
        pixel_object.save_as_picca_delta('tau',filename,header,notnorm=True,overwrite=args.overwrite,add_QSO_RSDs=args.add_QSO_RSDs,compress=args.compress,all_absorbers=args.picca_all_absorbers)

        #Picca flux
        filename = utils.get_file_name(location,'picca-flux-notnorm',N_side,pixel)
        pixel_object.save_as_picca_delta('flux',filename,header,notnorm=True,overwrite=args.overwrite,add_QSO_RSDs=args.add_QSO_RSDs,compress=args.compress,all_absorbers=args.picca_all_absorbers)

        """
        ## Disable this for the moment.
        #Save the final statistics file for this pixel.
        filename = utils.get_file_name(location,'statistics',N_side,pixel)
        statistics = pixel_object.save_statistics(filename,overwrite=args.overwrite,compress=args.compress,all_absorbers=args.picca_all_absorbers)
        """
    else:
        #If transmission_only is not False, remove the gaussian-colore file.
        os.remove(gaussian_filename)

    #print('{:3.2f} checkpoint RSD files'.format(time.time()-t)); t = time.time()

    return new_cosmology

#define the tasks
tasks = [(args.out_dir,pixel,args.nside,zero_mean_delta,args.lambda_min,args.tuning_file) for pixel in pixel_list]

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = args.nproc)
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
    header['NSIDE'] = args.nside

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
