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

parser.add_argument('--picca_N_merge', type = int, default = None, required=False,
                    help = 'number of cells to merge when rebinning for picca files')

args = parser.parse_args()

base_dir = args.base_dir
N_processes = args.nproc
N_side = args.nside
pixels = args.pixels
if not pixels:
    pixels = list(range(12*N_side**2))
N_merge = args.picca_N_merge

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

# TODO: paralellise
#Get the statistics from all pixels.
statistics_list = []
for pixel in pixels:
    s_filename = utils.get_file_name(base_dir,'statistics',N_side,pixel)
    s = fits.open(dirname+s_filename)
    statistics_list += s[1].data
    s.close()

#Combine the statistics from all of the pixels.
statistics = combine_statistics(statistics_list)

#Save the final file.
filename = 'statistics.fits'
stats.write_statistics(base_dir,filename,statistics)

################################################################################
"""
Renormalise the picca files using the global statistics data.
Also rebin the files to get the desired cell_size.
"""

#Also need to add in the rebinned ones? Or should we rebin here?
#Issue  with different names in the statistics files?
quantities = ['tau','flux']
stats_quantities = ['TAU','F']

# TODO: paralellise
#For each pixel, and each quantity, renormalise the picca file
for pixel in pixels:

    #Open up the per-pixel stats file
    dirname = utils.get_dir_name(base_dir,pixel)
    s_filename = utils.get_file_name(base_dir,'statistics',N_side,pixel)
    s = fits.open(dirname+s_filename)

    for i,q in enumerate(quantities):

        #Get the old mean, and renormalise.
        # TODO: this use of "stats_quantities" is v ugly
        old_mean = s[1].data[stats_quantities[i]+'_MEAN']
        filename = utils.get_file_name(base_dir,'picca-'+q,N_side,pixel)
        utils.renormalise_picca_file(filepath,old_mean,new_mean,N_merge=N_merge)
