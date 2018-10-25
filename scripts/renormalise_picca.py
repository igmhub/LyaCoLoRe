import numpy as np
from astropy.io import fits
from multiprocessing import Pool
import multiprocessing
import sys
import time

from pyacolore import tuning,utils

lya = utils.lya_rest

quantity = 'flux'
IVAR_cutoff = 1150.
basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/process_output_G_hZsmooth_4096_32_sr2.0_bm1_biasG18_picos_newNz_mpz0_seed1003_123_nside16/'
outdir = basedir
mean_data_filename = 'mean_data_{}_cut{}_v4.0.fits'.format(quantity,IVAR_cutoff)
min_number_cells = 2
rebin_size_hMpc = 2.3637
N_processes = 64

pixel_list = list(range(3072))

def renormalise(basedir,pixel,IVAR_cutoff,min_number_cells,rebin_size_hMpc,outdir):

    #Catalog dtype = [('RA', 'f8'), ('DEC', 'f8'), ('Z', 'f8'), ('Z_noRSD', 'f8'), ('MOCKID', int)]

    #Open the original picca flux file and extract the data.
    p_filename = basedir + '/{}/{}/picca-{}-16-{}.fits'.format(pixel//100,pixel,quantity,pixel)
    h = fits.open(p_filename)

    initial_delta = h[0].data.T
    IVAR = h[1].data.T
    LOGLAM_MAP = h[2].data
    catalog = h[3].data

    h.close()

    Z_QSO = catalog['Z']
    N_qso = initial_delta.shape[0]
    N_cells = initial_delta.shape[1]

    Z = (10**LOGLAM_MAP)/lya - 1
    
    #Load tha initial mean data.
    mean_data = fits.open(mean_data_filename)
    measured_mean_zs = mean_data[1].data['z']
    measured_mean = mean_data[1].data['mean']
    mean_data.close()

    #HACK
    measured_mean += 10 ** -10

    #Make new, proper delta rows.
    measured_mean = np.interp(Z,measured_mean_zs,measured_mean)
    new_delta = np.zeros(initial_delta.shape)
    for j in range(N_cells):
        new_delta[:,j] = (1 + initial_delta[:,j])/(1 + measured_mean[j]) - 1
    
    #new_delta = initial_delta
    new_delta *= IVAR

    #If needs be, rebin.
    if rebin_size_hMpc:

        #Open the gaussian colore file to get R.
        gc_filename = basedir + '/{}/{}/gaussian-colore-16-{}.fits'.format(pixel//100,pixel,pixel)
        gc = fits.open(gc_filename)
        Z_gc = gc[4].data['Z']
        R_gc = gc[4].data['R']
        gc.close()
        R = np.interp(Z,Z_gc,R_gc)

        #Rebin the relevant quantities
        grid_start = R[0]
        rebin_map = np.array((R - R[0])//rebin_size_hMpc,dtype='int')
        rebin_new_delta = np.zeros((N_qso,max(rebin_map)))
        rebin_IVAR = np.zeros((N_qso,max(rebin_map)))
        rebin_Z = np.zeros(max(rebin_map))
        j_hi = -1
        for j in range(max(rebin_map)):
            j_lo = j_hi + 1 #np.argmax(rebin_map==j)
            j_hi = np.searchsorted(rebin_map,j+1) - 1
            rebin_new_delta[:,j] = np.average(new_delta[:,j_lo:j_hi+1],axis=1)
            rebin_IVAR[:,j] = (np.sum(IVAR[:,j_lo:j_hi+1]) == j_hi+1 - j_lo)
            rebin_Z[j] = np.average(Z[j_lo:j_hi+1])
        rebin_LOGLAM_MAP = np.log10(lya*(1+rebin_Z))

        rebin_IVAR = utils.make_IVAR_rows(IVAR_cutoff,Z_QSO,rebin_LOGLAM_MAP)

        valid_cells = np.sum(rebin_IVAR,axis=1)
        relevant_QSOs = valid_cells > min_number_cells

        rebin_new_delta = rebin_new_delta[relevant_QSOs,:]
        rebin_IVAR = rebin_IVAR[relevant_QSOs,:]
        relevant_catalog = catalog[relevant_QSOs]

        picca_0 = rebin_new_delta.T
        picca_1 = rebin_IVAR.T
        picca_2 = rebin_LOGLAM_MAP
        picca_3 = relevant_catalog

        new_pf_filename = outdir + '/{}/{}/picca-{}-renorm-rebin-2.3637-16-{}.fits'.format(pixel//100,pixel,quantity,pixel)

    else:
        #Organise the data into picca-format arrays.
        picca_0 = new_delta.T
        picca_1 = IVAR.T
        picca_2 = LOGLAM_MAP
        picca_3 = catalog

        new_pf_filename = outdir + '/{}/{}/picca-{}-renorm-16-{}.fits'.format(pixel//100,pixel,quantity,pixel)

    header = h[0].header

    #Make the data into suitable HDUs.
    hdu_delta = fits.PrimaryHDU(data=picca_0,header=header)
    hdu_iv = fits.ImageHDU(data=picca_1,header=header,name='IV')
    hdu_LOGLAM_MAP = fits.ImageHDU(data=picca_2,header=header,name='LOGLAM_MAP')
    hdu_CATALOG = fits.BinTableHDU(data=picca_3,header=header,name='CATALOG')
    
    #Combine the HDUs into and HDUlist and save as a new file. Close the HDUlist.
    hdulist = fits.HDUList([hdu_delta, hdu_iv, hdu_LOGLAM_MAP, hdu_CATALOG])
    hdulist.writeto(new_pf_filename,overwrite=True)
    hdulist.close()

    return

################################################################################

"""
Define the multiprocessing tracking functions
"""

#Define a progress-tracking function.
def log_result(retval):

    results.append(retval)
    N_complete = len(results)
    N_tasks = len(tasks)

    utils.progress_bar(N_complete,N_tasks,start_time)

#Define an error-tracking function.
def log_error(retval):
    print('Error:',retval)

################################################################################

#define the tasks
tasks = [(basedir,pixel,IVAR_cutoff,min_number_cells,rebin_size_hMpc,outdir) for pixel in pixel_list]

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = N_processes)
    results = []
    start_time = time.time()

    for task in tasks:
        pool.apply_async(renormalise,task,callback=log_result,error_callback=log_error)

    pool.close()
    pool.join()

print('\nTime to renormalise files: {:4.0f}s.\n'.format(time.time()-start_time))
