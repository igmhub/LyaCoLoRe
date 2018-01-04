import numpy as np
from astropy.io import fits
from astropy.table import Table
import healpy as hp
import split_pixels_class
import time

#Set up a list of the files to import.
input_file_location = '/Users/jfarr/Projects/repixelise/test_input'
input_filename_structure = 'out_srcs_s0_{}.fits'
file_numbers = [16]
input_format = 'colore'

#Define the desired power of 2 for Nside for the output. This should be larger than that of the input.
N_side_pow2 = 2

N_side = 2**N_side_pow2
N_pix = 12*N_side**2

#Define the minimum value of z that we are interested in.
z_min = 1.85

#Define the directory that the output files should be saved to.
output_location = '/Users/jfarr/Projects/repixelise/test_output/WIP_class_output'
output_filename_structure = '{}_zmin_{}_nside_{}_pix_{}.fits'
output_format = 'picca'

#Define what the program should do if files with the proposed filename already exist. Options are:
#0 - save as a separate file with a suffix _v2 (and then _v3 if _v2 already exists etc.).
#1 - if there are duplicate quasars, follow #0. Otherwise merge the existing file with the new data.
#2 - Merge the existing file with the new data, avoiding duplicate objects.
existing_file_option = 0

start = time.time()
split_pixels_class.split_pixels(input_file_location,input_filename_structure,file_numbers,input_format,N_side,output_location,output_filename_structure,output_format,z_min=z_min,existing_file_option=existing_file_option)
end = time.time()
print('Class method takes: {}s'.format(end-start))
