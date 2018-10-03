import numpy as np
from astropy.io import fits
from multiprocessing import Pool
import multiprocessing
import sys
import time

from pyacolore import tuning,utils

lya = utils.lya_rest

IVAR_cutoff = 1150.
basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/process_output_G_hZsmooth_4096_32_sr2.0_bm1_biasG18_picos_newNz_mpz0_nside16'
min_number_cells = 2
rebin_size_hMpc = 3.5
N_processes = 32

pixel_list = list(range(3072))

def load_transmission_file(filename):

    h = fits.open(filename)
    CATALOG = h[1].data
    LAMBDA = h[2].data
    F = h[3].data
    h.close()

    return F, LAMBDA, CATALOG

def save_picca_flux_file():

    #Make the data into suitable HDUs.
    hdu_F = fits.PrimaryHDU(data=picca_0,header=header)
    hdu_iv = fits.ImageHDU(data=picca_1,header=header,name='IV')
    hdu_LOGLAM_MAP = fits.ImageHDU(data=picca_2,header=header,name='LOGLAM_MAP')
    cols_CATALOG = fits.ColDefs(picca_3)
    hdu_CATALOG = fits.BinTableHDU.from_columns(cols_CATALOG,header=header,name='CATALOG')

    #Combine the HDUs into and HDUlist and save as a new file. Close the HDUlist.
    hdulist = fits.HDUList([hdu_F, hdu_iv, hdu_LOGLAM_MAP, hdu_CATALOG])
    hdulist.writeto(filename)
    hdulist.close()

    return

def renorm_picca_flux(basedir,pixel,IVAR_cutoff,min_number_cells,rebin_size_hMpc):

    #F, LAMBDA, CATALOG = load_transmission_file(filename)
    #Catalog dtype = [('RA', 'f8'), ('DEC', 'f8'), ('Z', 'f8'), ('Z_noRSD', 'f8'), ('MOCKID', int)]

    #Open the original picca flux file and extract the data.
    pf_filename = basedir + '/{}/{}/picca-flux-16-{}.fits'.format(pixel//100,pixel,pixel)
    h = fits.open(pf_filename)

    initial_delta_F = h[0].data.T
    IVAR = h[1].data.T
    LOGLAM_MAP = h[2].data

    N_qso = initial_delta_F.shape[0]
    N_cells = initial_delta_F.shape[1]

    Z = (10**h[2].data)/lya - 1

    #Get the data for the mean_F originally used (i.e. from model).
    #Use it to go from delta_F to F.
    initial_mean_F = tuning.get_mean_F_model(Z)
    F = np.zeros(initial_delta_F.shape)
    for j in range(N_cells):
        F[:,j] = (initial_delta_F[:,j] + 1) * initial_mean_F[j]

    #Load tha actual mean_F data.
    mean_data = fits.open('mean_data.fits')
    new_mean_F_zs = mean_data[1].data['z']
    new_mean_F = mean_data[1].data['mean']
    mean_data.close()

    #HACK
    new_mean_F += 10 ** -10

    #Make new, proper delta_F rows.
    new_mean_F = np.interp(Z,new_mean_F_zs,new_mean_F)
    new_delta_F = np.zeros(initial_delta_F.shape)
    for j in range(N_cells):
        new_delta_F[:,j] = (F[:,j]/new_mean_F[j]) - 1

    new_delta_F *= IVAR

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
        rebin_new_delta_F = np.zeros((N_qso,max(rebin_map)))
        rebin_IVAR = np.zeros((N_qso,max(rebin_map)))
        rebin_Z = np.zeros(max(rebin_map))
        j_hi = -1
        for j in range(max(rebin_map)):
            j_lo = j_hi + 1 #np.argmax(rebin_map==j)
            j_hi = np.searchsorted(rebin_map,j+1) - 1
            rebin_new_delta_F[:,j] = np.average(new_delta_F[:,j_lo:j_hi+1],axis=1)
            rebin_IVAR[:,j] = (np.sum(IVAR[:,j_lo:j_hi+1]) == j_hi+1 - j_lo)
            rebin_Z[j] = np.average(Z[j_lo:j_hi+1])        
        rebin_LOGLAM_MAP = np.log10(lya*(1+rebin_Z))

        picca_0 = rebin_new_delta_F.T
        picca_1 = rebin_IVAR.T
        picca_2 = rebin_LOGLAM_MAP

        new_pf_filename = basedir + '/{}/{}/picca-flux-renorm-rebin-16-{}.fits'.format(pixel//100,pixel,pixel)

    else:
        #Organise the data into picca-format arrays.
        picca_0 = new_delta_F.T
        picca_1 = IVAR.T
        picca_2 = LOGLAM_MAP

        new_pf_filename = basedir + '/{}/{}/picca-flux-renorm-16-{}.fits'.format(pixel//100,pixel,pixel)

    picca_3 = h[3].data
    header = h[0].header

    #Make the data into suitable HDUs.
    hdu_F = fits.PrimaryHDU(data=picca_0,header=header)
    hdu_iv = fits.ImageHDU(data=picca_1,header=header,name='IV')
    hdu_LOGLAM_MAP = fits.ImageHDU(data=picca_2,header=header,name='LOGLAM_MAP')
    cols_CATALOG = fits.ColDefs(picca_3)
    hdu_CATALOG = fits.BinTableHDU.from_columns(cols_CATALOG,header=header,name='CATALOG')

    #Combine the HDUs into and HDUlist and save as a new file. Close the HDUlist.
    hdulist = fits.HDUList([hdu_F, hdu_iv, hdu_LOGLAM_MAP, hdu_CATALOG])
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
tasks = [(basedir,pixel,IVAR_cutoff,min_number_cells,rebin_size_hMpc) for pixel in pixel_list]

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = N_processes)
    results = []
    start_time = time.time()

    for task in tasks:
        pool.apply_async(renorm_picca_flux,task,callback=log_result,error_callback=log_error)

    pool.close()
    pool.join()

print('\nTime to renormalise files: {:4.0f}s.\n'.format(time.time()-start_time))
