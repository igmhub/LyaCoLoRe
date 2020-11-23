#!/usr/bin/env python

import multiprocessing
import numpy as np
import os
import sys
import time

from astropy.io import fits
from multiprocessing import Pool
from scipy.interpolate import interp1d

from lyacolore import catalog, parse, simulation_data, transformation, utils

################################################################################

#Script to produce a a set of transmission files from CoLoRe's output files.

# TODO: Tidy up measuring of SIGMA_G and subsequent DLA method.
# TODO: Exchange lambda_min for z_min for cells.

################################################################################

"""
Set up the file locations and filename structures.
Also define option preferences.
"""

args = parse.get_args(sys.argv)

################################################################################
"""
Check that the arguments are all ok, and construct the MOCKID_lookup from the
master file.
"""

#Get the simulation parameters from the parameter file.
simulation_parameters = utils.get_simulation_parameters(args.in_dir,args.param_file)

with fits.open(args.out_dir+'/master.fits') as master:
    MOCKID_lookup = catalog.make_MOCKID_lookup(master[1].data,args.pixels)
pixel_list = list(sorted(set([k[1] for k in MOCKID_lookup.keys()])))

################################################################################
"""
Construct pixel objects.
"""

print('Working on per-HEALPix pixel initial skewer files...')
pixel_objects = simulation_data.get_pixel_objects(pixel_list,args,MOCKID_lookup)

################################################################################
"""
To correctly calculate the physical fields when using Gaussian input skewers, we
must combine the values from each pixel, to compute an overall value.
"""

if args.skewer_type == 'gaussian':
    w = np.sum([p.IVAR_rows.astype(int).sum() for p in pixel_objects.values()])
    wsg2 = np.sum([p.IVAR_rows.astype(int).sum()*(p.SIGMA_G**2) for p in pixel_objects.values()])
    sigma_g_global = np.sqrt(wsg2/w)

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

def produce_final_skewers(pixel,args):

    ## Extract the pixel object.
    pixel_object = pixel_objects[pixel]

    ## Modify the value of sigma_g appropriately.
    if args.skewer_type == 'gaussian':
        pixel_object.SIGMA_G = sigma_g_global

    simulation_data.process_skewers(pixel_object,pixel,args,save_stages=(~args.transmission_only),save_transmission=True)

    #Make a variable containing the new cosmology data.
    new_cosmology = pixel_object.return_cosmology()

    return new_cosmology

#define the tasks
tasks = [(pixel,args) for pixel in pixel_list]

#Run the multiprocessing pool
if __name__ == '__main__':
    results = utils.run_multiprocessing(produce_final_skewers,tasks,args.nproc)

################################################################################
"""
PROBABLY COULD MOVE THIS TO make_summaries
Having added small scale power, we must add a new HDU to the master file's cosmology.
"""

print('\nUpdating master file\'s cosmology...')
#First check that the new cosmologies are all the same.
# TODO: some kind of system to check consistency here?
new_cosmology = results[0]

#Reorganise the data.
with fits.open(args.out_dir+'/master.fits',mode='update') as m:
    if len(master)<=3:

        ## Make clear that this is the CoLoRe cosmology.
        m['COSMO'].name = 'COSMO_COL'

        ## Make an appropriate header.
        header = fits.Header()
        header['NSIDE'] = args.nside

        ## Add the expanded cosmology.
        hdu_cosmo_exp = fits.BinTableHDU.from_columns(new_cosmology,header=header,name='COSMO_EXP')
        m.append(hdu_cosmo_exp)

        ## Update the file.
        m.flush()

print('Process complete!\n')

################################################################################

"""
Celebrate!
"""
