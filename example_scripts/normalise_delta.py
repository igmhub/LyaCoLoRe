import numpy as np
from astropy.io import fits

#Script to normalise the value of delta such that its mean is 0.

#Set up a list of the files to import.
file_location = '/Users/jfarr/Projects/repixelise/test_input'
file_prefix = 'out_srcs_s0_'
file_numbers = [16]
files = split.file_list(file_location,file_prefix,file_numbers)

for i, file in enumerate files:
    initial = fits.open(file)
    DELTA = initial[2].data
    new_DENSITY = np.array([1+x for x in DELTA])/(1+np.mean(DELTA))
    new_DELTA = np.array([x-1 for x in new_DENSITY])

    final = initial
    final[2].data = new_DELTA

    new_filename = file_location + file_prefix + 'normalised_delta_' + file_numbers[i] + '.fits'
    final.write(new_filename)
