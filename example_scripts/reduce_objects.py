import numpy as np
from astropy.io import fits
import cut_fits as cut

#Reduces the number of objects in a CoLoRe output file.

#Input the original filename and location.
original_filename='/Users/jfarr/Projects/repixelise/test_input/out_srcs_s0_0.fits'

#Set the number of objects desired in the output.
N_obj_desired = 5000

#Specify the mode of object selection: 0=select the last (highest z) N_obj_desired objects, 1=select N_obj_desired objects at random
mode = 1

#Cut the file
cut.cut_fits(original_filename, N_obj_desired, mode)
