import numpy as np
from astropy.io import fits

#Function to normalise the value of delta from a CoLoRe output file such that its mean is 0.

def normalise_delta(file_location,file_prefix,file_number):

    file = file_location + '/' + file_prefix + str(file_numbers[i]) + '.fits'
    initial = fits.open(file)
    DELTA = initial[2].data
    new_DENSITY = np.array([1+x for x in DELTA])/(1+np.mean(DELTA))
    new_DELTA = np.array([x-1 for x in new_DENSITY])

    final = initial
    final[2].data = new_DELTA

    new_filename = file_location + '/' + file_prefix + 'normalised_delta_' + str(file_numbers[i]) + '.fits'
    final.writeto(new_filename)

    initial.close()

    return
