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

#Make the DLA master file
DLA.make_DLA_master(base_dir,N_side,pixels)

################################################################################
"""
Make the global statistics file.
"""

print('Making the global statistics file...')

# TODO: paralellise
#Get the statistics from all pixels.
statistics_list = []

def get_statistics(pixel):
    dirname = utils.get_dir_name(base_dir,pixel)

    #Open up the statistics file without RSDs and extract data.
    s_noRSD_filename = utils.get_file_name(dirname,'statistics-noRSD',N_side,pixel)
    s_noRSD = fits.open(s_noRSD_filename)
    statistics_noRSD = s_noRSD[1].data
    s_noRSD.close()

    #Open up the statistics file with RSDs and extract data.
    s_filename = utils.get_file_name(dirname,'statistics',N_side,pixel)
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
stats.write_statistics(base_dir,filename,statistics_noRSD)

statistics = stats.combine_statistics(statistics_list)
filename = './statistics.fits'
stats.write_statistics(base_dir,filename,statistics)

################################################################################
"""
We need to renormalise the picca files using global statistics data.
First we create a global file and save it.
We then ensure that the global mean of all the picca-deltas is 0.
Also, we rebin the files to get the desired cell_size.
"""

print('Renormalising the picca files...')

#Also need to add in the rebinned ones? Or should we rebin here?
#Issue  with different names in the statistics files?
quantities = ['tau','flux']
stats_quantities = ['TAU','F']

tasks = [(pixel,N_merge) for pixel in pixels for N_merge in N_merge_values]

#For each pixel, and each quantity, renormalise the picca file
def renormalise(pixel,N_merge):
    #Open up the per-pixel stats file
    dirname = utils.get_dir_name(base_dir,pixel)
    s_filename = utils.get_file_name(dirname,'statistics',N_side,pixel)
    s = fits.open(s_filename)

    for i,q in enumerate(quantities):

        #Get the old mean, and renormalise.
        # TODO: this use of "stats_quantities" is v ugly
        lookup_name = stats_quantities[i]+'_MEAN'
        #old_mean = s[1].data[lookup_name]
        old_mean = np.ones(s[1].data.shape)

        #Renormalise the files without RSDs.
        new_mean = statistics_noRSD[lookup_name]
        filename = utils.get_file_name(dirname,'picca-'+q+'-noRSD',N_side,pixel)
        if N_merge == 1:
            out = utils.get_file_name(dirname,'picca-'+q+'-noRSD-renorm-',N_side,pixel)
        else:
            out = utils.get_file_name(dirname,'picca-'+q+'-noRSD-renorm-rebin-{}-'.format(N_merge)+q,N_side,pixel)
        utils.renormalise_picca_file(filename,old_mean,new_mean,N_merge=N_merge,out_filepath=out)

        #Renormalise the files with RSDs.
        new_mean = statistics[lookup_name]
        filename = utils.get_file_name(dirname,'picca-'+q,N_side,pixel)
        if N_merge == 1:
            out = utils.get_file_name(dirname,'picca-'+q+'-renorm-',N_side,pixel)
        else:
            out = utils.get_file_name(dirname,'picca-'+q+'-renorm-rebin-{}-'.format(N_merge)+q,N_side,pixel)
        utils.renormalise_picca_file(filename,old_mean,new_mean,N_merge=N_merge,out_filepath=out)

    return

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = N_processes)
    results = []
    start_time = time.time()

    for task in tasks:
        pool.apply_async(renormalise,task,callback=log_result,error_callback=log_error)

    pool.close()
    pool.join()

################################################################################

"""
Celebrate!
"""
