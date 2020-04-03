#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from multiprocessing import Pool
import multiprocessing
import sys
import time
import os
import glob
import argparse

from lyacolore import utils, independent, stats, convert, simulation_data, DLA, RSD, tuning

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--base-dir', type = str, default = None, required=True,
                    help = 'directory of LyaCoLoRe output')

parser.add_argument('--nproc', type = int, default = 1, required=False,
                    help = 'number of processes to use')

parser.add_argument('--nside', type = int, default = 16, required=False,
                    help = 'HEALPix nside for output files (must be 2^n)')

parser.add_argument('--pixels', type = int, default = None, required=False,
                    help = 'which pixel numbers to work on', nargs='*')

parser.add_argument('--picca-N-merge-values', type = int, default = [1], required = False,
                    help = 'number of cells to merge when rebinning for picca files', nargs='*')

parser.add_argument('--overwrite', action="store_true", default = False, required=False,
                    help = 'overwrite existing files')

parser.add_argument('--compressed-input', action="store_true", default = False, required=False,
                    help = 'compress output files to .fits.gz')

parser.add_argument('--compress', action="store_true", default = False, required=False,
                    help = 'compress output files to .fits.gz')

parser.add_argument('--transmission-only', action="store_true", default = False, required=False,
                    help = 'save only the transmission file')

parser.add_argument('--make-DLA-master', action="store_true", default = False, required=False,
                    help = 'make a DLA master file')

parser.add_argument('--add-picca-DLA-drqs', action="store_true", default = False, required=False,
                    help = 'save picca format drq files for DLAs')

parser.add_argument('--footprint', type = str, default = None, required = False,
                    choices=['full_sky','desi','desi_pixel','desi_pixel_plus'],
                    help = 'name of footprint to use')

args = parser.parse_args()

base_dir = args.base_dir
N_processes = args.nproc
N_side = args.nside
pixels = args.pixels
if not pixels:
    pixels = list(range(12*N_side**2))
N_merge_values = args.picca_N_merge_values
if not N_merge_values:
    N_merge_values = []
overwrite = args.overwrite
compressed_input = args.compressed_input
compress = args.compress
transmission_only = args.transmission_only
add_picca_DLA_drqs = args.add_picca_DLA_drqs

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
Make the DLA master file.
"""

if args.make_DLA_master:
    print('Making the DLA master file...')

    def get_DLA_data(pixel):
        dirname = utils.get_dir_name(base_dir,pixel)
        filename = utils.get_file_name(dirname,'transmission',N_side,pixel,compressed=compressed_input)
        DLA_data = DLA.get_DLA_data_from_transmission(pixel,filename)
        return DLA_data

    if args.footprint in ['desi','desi_pixel']:
        footprint_pixels = np.loadtxt('/global/homes/j/jfarr/Projects/LyaCoLoRe/input_files/DESI_pixels.txt')
    elif args.footprint in ['desi_pixel_plus']:
        footprint_pixels = np.loadtxt('/global/homes/j/jfarr/Projects/LyaCoLoRe/input_files/DESI_pixels_plus.txt')
    else:
        footprint_pixels = np.array(list(range(12*args.nside*args.nside)))

    tasks = [(pixel,) for pixel in pixels if pixel in footprint_pixels]

    #Run the multiprocessing pool
    if __name__ == '__main__':
        pool = Pool(processes = N_processes)
        results = []
        start_time = time.time()

        for task in tasks:
            pool.apply_async(get_DLA_data,task,callback=log_result,error_callback=log_error)

        pool.close()
        pool.join()

    #Make the DLA master file
    filename = base_dir + '/master_DLA.fits'
    DLA.write_DLA_master(results,filename,N_side,overwrite=overwrite)

    #Write DLA data to picca DRQ format files if desired.
    if add_picca_DLA_drqs:
        DLA.write_DLA_master(results,base_dir+'/master_DLA_picca_RSD.fits',N_side,overwrite=overwrite,picca_DRQ=True,add_DRQ_RSDs=True)
        DLA.write_DLA_master(results,base_dir+'/master_DLA_picca_NO_RSD.fits',N_side,overwrite=overwrite,picca_DRQ=True,add_DRQ_RSDs=False)

################################################################################
"""
Make the global statistics file.
"""

if not transmission_only:
    print('Making the global statistics file...')

# TODO: paralellise
    #Get the statistics from all pixels.
    statistics_list = []

    def get_statistics(pixel):
        dirname = utils.get_dir_name(base_dir,pixel)

        #Open up the statistics file without RSDs and extract data.
        s_noRSD_filename = utils.get_file_name(dirname,'statistics-noRSD',N_side,pixel,compressed=compressed_input)
        s_noRSD = fits.open(s_noRSD_filename)
        statistics_noRSD = s_noRSD[1].data
        s_noRSD.close()

        #Open up the statistics file with RSDs and extract data.
        s_filename = utils.get_file_name(dirname,'statistics',N_side,pixel,compressed=compressed_input)
        s = fits.open(s_filename)
        statistics = s[1].data
        s.close()

        return statistics_noRSD, statistics

    tasks = [(pixel,) for pixel in pixels]

    #Run the multiprocessing pool
    if __name__ == '__main__':
        pool = Pool(processes = N_processes)
        results = []
        start_time = time.time()

        for task in tasks:
            pool.apply_async(get_statistics,task,callback=log_result,error_callback=log_error)

        pool.close()
        pool.join()

    #Make lists with and without RSDs
    statistics_noRSD_list = []
    statistics_list = []
    for result in results:
        statistics_noRSD_list += [result[0]]
        statistics_list += [result[1]]

    #Combine the statistics from all of the pixels and save, with and without RSDs.
    statistics_noRSD = stats.combine_statistics(statistics_noRSD_list)
    filename = './statistics_noRSD.fits'
    stats.write_statistics(base_dir+filename,statistics_noRSD,overwrite=overwrite)

    statistics = stats.combine_statistics(statistics_list)
    filename = './statistics.fits'
    stats.write_statistics(base_dir+filename,statistics,overwrite=overwrite)

################################################################################
"""
We need to renormalise the picca files using global statistics data.
First we create a global file and save it.
We then ensure that the global mean of all the picca-deltas is 0.
Also, we rebin the files to get the desired cell_size.
"""

if not transmission_only:
    print('Renormalising and rebinning the picca files...')

    #Also need to add in the rebinned ones? Or should we rebin here?
    #Issue  with different names in the statistics files?
    type_1_quantities = ['gaussian','density']
    type_2_quantities = ['tau','flux']
    stats_quantities = ['TAU','F']

    tasks = [(pixel,) for pixel in pixels]

    #For each pixel, and each quantity, renormalise the picca file
    def normalise_and_rebin(pixel):

        #Get the directory name.
        dirname = utils.get_dir_name(base_dir,pixel)

        for N_merge in N_merge_values:
            for i,q in enumerate(type_1_quantities):

                #print('rebinning {} file'.format(q))
                t = time.time()
                #Rebin the files.
                filename = utils.get_file_name(dirname,'picca-'+q,N_side,pixel,compressed=compressed_input)
                if N_merge > 1:
                    out = utils.get_file_name(dirname,'picca-'+q+'-rebin-{}'.format(N_merge),N_side,pixel)
                    utils.renorm_rebin_picca_file(filename,N_merge=N_merge,out_filepath=out,overwrite=overwrite)
                #print('--> {:1.3f}s'.format(time.time()-t))

            for i,q in enumerate(type_2_quantities):

                #print('getting stats data')
                t = time.time()
                #Get the old mean, and renormalise.
                # TODO: this use of "stats_quantities" is v ugly
                lookup_name = stats_quantities[i]+'_MEAN'
                #print('--> {:1.3f}s'.format(time.time()-t))

                #print('rebin/renorm-ing {} noRSD file'.format(q))
                t = time.time()
                #Renormalise the files without RSDs.
                #old_mean = s_noRSD[1].data[lookup_name]
                old_mean = None
                new_mean = statistics_noRSD[lookup_name]
                filename = utils.get_file_name(dirname,'picca-'+q+'-noRSD-notnorm',N_side,pixel,compressed=compressed_input)
                if N_merge == 1:
                    out = utils.get_file_name(dirname,'picca-'+q+'-noRSD',N_side,pixel)
                else:
                    out = utils.get_file_name(dirname,'picca-'+q+'-noRSD-rebin-{}'.format(N_merge),N_side,pixel)
                #print(out)
                utils.renorm_rebin_picca_file(filename,old_mean=old_mean,new_mean=new_mean,N_merge=N_merge,out_filepath=out,overwrite=overwrite,compress=compress)
                #print('--> {:1.3f}s'.format(time.time()-t))

                #print('rebin/renorm-ing {} RSD file'.format(q))
                t = time.time()
                #Renormalise the files with RSDs.
                #old_mean = s[1].data[lookup_name]
                old_mean = None
                new_mean = statistics[lookup_name]
                filename = utils.get_file_name(dirname,'picca-'+q+'-notnorm',N_side,pixel,compressed=compressed_input)
                if N_merge == 1:
                    out = utils.get_file_name(dirname,'picca-'+q,N_side,pixel)
                else:
                    out = utils.get_file_name(dirname,'picca-'+q+'-rebin-{}'.format(N_merge),N_side,pixel)
                #print(out)
                utils.renorm_rebin_picca_file(filename,old_mean=old_mean,new_mean=new_mean,N_merge=N_merge,out_filepath=out,overwrite=overwrite,compress=compress)
                #print('--> {:1.3f}s'.format(time.time()-t))

        return

    #Run the multiprocessing pool
    if __name__ == '__main__':
        pool = Pool(processes = N_processes)
        results = []
        start_time = time.time()

        for task in tasks:
            pool.apply_async(normalise_and_rebin,task,callback=log_result,error_callback=log_error)

        pool.close()
        pool.join()

################################################################################
"""
Compress all output if desired.
"""

if compress:
    print('Compressing output files...')

    compress_file = utils.compress_file
    fi = glob.glob(base_dir+'/*/*/*.fits')
    tasks = [(f,) for f in fi]

    #Run the multiprocessing pool
    if __name__ == '__main__':
        pool = Pool(processes = N_processes)
        results = []
        start_time = time.time()

        for task in tasks:
            pool.apply_async(compress_file,task,callback=log_result,error_callback=log_error)

        pool.close()
        pool.join()

################################################################################

"""
Celebrate!
"""
