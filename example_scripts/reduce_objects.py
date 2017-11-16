import numpy as np
from astropy.io import fits
import cut_fits as cut

#Reduces the number of objects in a CoLoRe output file.

#Input the original filename and location.
file_location = '/Users/jfarr/Projects/repixelise/test_input'
file_name = 'out_srcs_s0_0.fits'

#Set the number of objects desired in the output.
N_obj_desired = 1000

#Specify the mode of object selection: 0=select the last (highest z) N_obj_desired objects, 1=select N_obj_desired objects at random
mode = 1

#Cut the file
cut.cut_fits(file_location, file_name, N_obj_desired, mode)
