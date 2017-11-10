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
#ifiles = [file_location+filename_1, file_location+filename_2]

#Option 2: set up a prefix/suffix structure.
#Set up a list of the files to import.
file_location = '/Users/jfarr/Projects/repixelise/test_input'
file_prefix = 'test_skewers_4096_gaussian_srcs_s0_'
file_numbers = [0]
files = split.file_list(file_location,file_prefix,file_numbers)

#Define the desired power of 2 for Nside for the output. This should be larger than that of the input.
Nside_pow2 = 1
Nside = 2**Nside_pow2
Npix = 12*Nside**2

#Define the directory that the output files should be saved to.
save_location = '/Users/jfarr/Projects/repixelise/test_output/'

#Split the files.
split.split_pixels(Nside,files,save_location)
