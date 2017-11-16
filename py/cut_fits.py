import numpy as np
from astropy.io import fits
import random

# Takes a CoLoRe output file, and cuts it to have N_obj objects.
def cut_fits(original_filename, N_obj_desired, mode):

    original = fits.open(original_filename)
    N_obj_original = original[1].data.shape[0]

    if N_obj_desired > N_obj_original:
        print('There are only %d objects in the file, so you cannot select %d of them.' % N_obj_original, N_obj_desired)
        quit()
    else:
        if mode == 0:

            # Option 0: select the last N_obj_desired entries.
            new_1 = original[1].data[-N_obj_desired:]
            new_2 = original[2].data[-N_obj_desired:]
            new_3 = original[3].data[-N_obj_desired:]

        elif mode == 1:

            # Option 1: select a random set of N_obj_desired entries.
            sample = sorted(random.sample(range(N_obj_original),N_obj_desired))

            new_1 = original[1].data[sample]
            new_2 = original[2].data[sample,:]
            new_3 = original[3].data[sample,:]
            new_4 = original[4].data


        else:
            print('Mode %d not recognised; please try again.' % mode)


        prihdu = original[0]
        hdu_new_1 = fits.BinTableHDU(new_1,name='Catalog')
        hdu_new_2 = fits.ImageHDU(data=new_2,header=None,name='Density skewers')
        hdu_new_3 = fits.ImageHDU(data=new_3,header=None,name='Velocity skewers')
        hdu_new_4 = original[4]

        hdulist = fits.HDUList([prihdu, hdu_new_1, hdu_new_2, hdu_new_3, hdu_new_4])

        #Re write this to give a better filename (use prefix and suffix)
        original_filename_nofits = original_filename[:-5]
        hdulist.writeto(original_filename_nofits + 'cut_%d' % N_obj_desired + '.fits')

        return
