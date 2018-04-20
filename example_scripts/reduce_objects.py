import numpy as np
from astropy.io import fits
import cut_fits as cut
import sys

#Reduces the number of objects in a CoLoRe output file.

#Input the original filename and location.
file_location = '/Users/James/Projects/test_data/output_G_hZ_4096_32_sr2.0_bm1/'
file_prefix = 'out_srcs_s1_'
file_suffices = [0]

#Specify the structure of the original file:
#0 = CoLoRe output format
#1 = picca input format
structure = 0

#Set the number of objects desired in the output.
N_obj_desired = int(sys.argv[1])

#Specify the mode of object selection:
#0=select the last (highest z) N_obj_desired objects
#1=select N_obj_desired objects at random
mode = 1

#Cut the file
for suffix in file_suffices:
    print('Working on {}.'.format(file_prefix+str(suffix)+'.fits'))
    cut.cut_fits(file_location, file_prefix+str(suffix)+'.fits', structure, N_obj_desired, mode)
