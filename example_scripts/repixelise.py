import numpy as np
from astropy.io import fits
from astropy.table import Table
import healpy as hp
import split_pixels as split

#Define the files that are to be repixelised.
#Option 1: name each file individually.
#file_location =
#filename_1 =
#filename_2 =
#files = [file_location+filename_1, file_location+filename_2]

#Option 2: set up a prefix/suffix structure.
#Set up a list of the files to import.
file_location = '/Users/jfarr/Projects/repixelise/test_input'
file_prefix = 'out_srcs_s0_'
file_numbers = [15]
files = split.file_list(file_location,file_prefix,file_numbers)

#Define the desired power of 2 for Nside for the output. This should be larger than that of the input.
N_side_pow2 = 2
N_side = 2**N_side_pow2
N_pix = 12*N_side**2

#Determine which format is desired for output. 0=from_fitsio, 1=from_image
#NOTE: Only output in the 'from_image' form is working now.
output_format = 1

#Define the directory that the output files should be saved to.
save_location = '/Users/jfarr/Projects/repixelise/test_output'

#Split the files.
split.split_pixels(N_side,files,file_numbers,save_location,output_format)
