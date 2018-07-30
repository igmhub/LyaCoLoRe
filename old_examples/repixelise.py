import numpy as np
from astropy.io import fits
from astropy.table import Table
import healpy as hp
import split_pixels as split
import WIP_split_pixels as WIP
import WIP_split_pixels_class as WIP_class
import time

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
input_format = 'colore'
file_numbers = [0]

#Define the desired power of 2 for Nside for the output. This should be larger than that of the input.
N_side_pow2 = 2

N_side = 2**N_side_pow2
N_pix = 12*N_side**2

#Define the minimum value of z that we are interested in.
z_min = 1.85

#Define the directory that the output files should be saved to.
save_location = '/Users/jfarr/Projects/repixelise/test_output'
new_filename_structure = 'zmin_{}_nside_{}_pix_{}.fits'
output_format = 'picca'

start = time.time()
WIP.split_pixels(original_file_location,original_filename_structure,file_numbers,input_format,N_side,z_min,save_location + '/WIP_output',new_filename_structure,output_format)
end = time.time()
print("WIP method: ",end-start)

#start = time.time()
#WIP_class.split_pixels(original_file_location,original_filename_structure,file_numbers,N_side,z_min,save_location+'/WIP_class_output',new_filename_structure,output_format)
#end = time.time()
#print("WIP class method: ",end-start)
