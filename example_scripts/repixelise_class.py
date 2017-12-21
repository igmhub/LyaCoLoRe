import numpy as np
from astropy.io import fits
from astropy.table import Table
import healpy as hp
import WIP_split_pixels_class as WIP_class
import time

"""
split_pixels(N_side,original_location,file_numbers,save_location,output_format,z_min)
"""

#Define the files that are to be repixelised.
#Option 1: name each file individually.
#file_location =
#filename_1 =
#filename_2 =
#files = [file_location+filename_1, file_location+filename_2]

#Option 2: set up a prefix/suffix structure.
#Set up a list of the files to import.
original_file_location = '/Users/jfarr/Projects/repixelise/test_input'
original_filename_structure = 'out_srcs_s0_{}.fits'
file_numbers = [16]
input_format = 'colore'

#Define the desired power of 2 for Nside for the output. This should be larger than that of the input.
N_side_pow2 = 2

N_side = 2**N_side_pow2
N_pix = 12*N_side**2

#Define the minimum value of z that we are interested in.
z_min = 1.85

#Define the directory that the output files should be saved to.
save_location = '/Users/jfarr/Projects/repixelise/test_output/WIP_class_output'
new_filename_structure = '{}_zmin_{}_nside_{}_pix_{}.fits'
output_format = 'picca'

#Define what the program should do if files with the proposed filename already exist. Options are:
#0 - save as a separate file with a suffix _2 (and then _3 if _2 already exists etc.).
#1 - if there are duplicate quasars, follow #0. Otherwise merge the existing file with the new data.
#2 - Merge the existing file with the new data, avoiding duplicate objects.
existing_file_option = 2

start = time.time()
WIP_class.split_pixels(original_file_location,original_filename_structure,file_numbers,input_format,N_side,save_location,new_filename_structure,output_format,z_min=z_min,existing_file_option=existing_file_option)
end = time.time()
print('Class method takes: {}s'.format(end-start))
