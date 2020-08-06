import numpy as np
from astropy.io import fits
import random

# Takes a CoLoRe output file, and cuts it to have N_obj objects.
def cut_fits(file_location, file_name, structure, N_obj_desired, mode):

    original_filename = file_location + '/' + file_name
    original = fits.open(original_filename)

    if structure == 0:
        N_obj_original = original[1].data.shape[0]

        if N_obj_desired > N_obj_original:
            print('There are only {} objects in the file, so you cannot select {} of them.'.format(N_obj_original, N_obj_desired))
            quit()

        else:
            if mode == 0:

                new_1 = original[1].data[-N_obj_desired:]
                new_2 = original[2].data[-N_obj_desired:]
                new_3 = original[3].data[-N_obj_desired:]

            elif mode == 1:

                sample = sorted(random.sample(range(N_obj_original),N_obj_desired))

                new_1 = original[1].data[sample]
                new_2 = original[2].data[sample,:]
                new_3 = original[3].data[sample,:]
                new_4 = original[4].data

                headers_1 = original[1].header
                headers_2 = original[2].header
                headers_3 = original[3].header
                headers_4 = original[4].header

                #Check if there's now irrelevant cells (i.e. all zeros). If so, trim them.
                #sum_of_new_2_columns = np.sum(new_2,axis=0)
                #zero_sum_start = np.argmax(sum_of_new_2_columns==0)

                #if zero_sum_start==0 and sum_of_new_2_columns[0]==0:
                #    exit()
                #elif zero_sum_start>0:
                #    new_2 = new_2[:,0:zero_sum_start]
                #    new_3 = new_3[:,0:zero_sum_start]
                #    new_4 = new_4[0:zero_sum_start]

            else:
                print('Mode {} not recognised; please try again.'.format(mode))


            prihdu = original[0]
            hdu_new_1 = fits.BinTableHDU(new_1,header=headers_1,name='Catalog')
            hdu_new_2 = fits.ImageHDU(data=new_2,header=headers_2,name='Density skewers')
            hdu_new_3 = fits.ImageHDU(data=new_3,header=headers_3,name='Velocity skewers')
            hdu_new_4 = fits.BinTableHDU(new_4,header=headers_4,name='Background cosmology')

            hdulist = fits.HDUList([prihdu, hdu_new_1, hdu_new_2, hdu_new_3, hdu_new_4])

            new_filename = file_location + '/' + 'N{}_'.format(N_obj_desired) + file_name

            hdulist.writeto(new_filename)


    elif structure == 1:
        N_obj_original = original[0].data.shape[1]

        if N_obj_desired > N_obj_original:
            print('There are only {} objects in the file, so you cannot select {} of them.'.format(N_obj_original, N_obj_desired))
            quit()

        else:
            if mode == 0:

                # Option 0: select the last N_obj_desired entries.
                new_0 = original[0].data[:,-N_obj_desired:]
                new_1 = original[1].data[:,-N_obj_desired:]
                new_3 = original[3].data[-N_obj_desired:]

            elif mode == 1:

                # Option 1: select a random set of N_obj_desired entries.
                # Sometimes these seem to hang, but I'm not sure why.
                sample = sorted(random.sample(range(N_obj_original),N_obj_desired))
                new_0 = original[0].data[:,sample]
                new_1 = original[1].data[:,sample]
                new_2 = original[2].data
                new_3 = original[3].data[sample]


                #Check if there's now irrelevant cells (i.e. all zeros). If so, trim them.
                sum_of_new_0_rows = np.sum(new_0,axis=1)
                zero_sum_start = np.argmax(sum_of_new_0_rows==0)

                if zero_sum_start==0 and sum_of_new_0_rows[0]==0:
                    exit()
                elif zero_sum_start>0:
                    new_0 = new_0[0:zero_sum_start,:]
                    new_1 = new_1[0:zero_sum_start,:]
                    new_2 = new_2[0:zero_sum_start]

            else:
                print('Mode {} not recognised; please try again.'.format(mode))

            prihdu = fits.PrimaryHDU(new_0)
            hdu_new_1 = fits.ImageHDU(data=new_1,header=None,name='IV')
            hdu_new_2 = fits.ImageHDU(data=new_2,header=None,name='LOGLAM_MAP')
            hdu_new_3 = fits.BinTableHDU(new_3,name='CATALOG')

            hdulist = fits.HDUList([prihdu, hdu_new_1, hdu_new_2, hdu_new_3])

            new_filename = file_location + '/' + 'N{}_'.format(N_obj_desired) + file_name

            hdulist.writeto(new_filename)

    else:
        print('Structure {} not recognised, please try again.'.format(structure))

    original.close()

    return
