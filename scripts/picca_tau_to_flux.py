import numpy as np
import glob
import time
from multiprocessing import Pool
import multiprocessing
from astropy.io import fits

basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v5/test_runs_with_bias/test_a2.0_b1.65_wb/'
N_merge = 10
min_number_cells = 2
pixels = np.array(list(range(200)))

statistics_file = basedir + '/statistics.fits'

def convert_picca_tau_to_flux(pixel):
    #for each pixel
    dirname = utils.get_dir_name(basedir,pixel)
    filepath = utils.get_file_name(dirname,'picca-tau-notnorm',16,pixel)
    h = fits.open(filepath)

    #Get data
    skewer_rows = h[0].data.T
    IVAR_rows = h[1].data.T
    LOGLAM_MAP = h[2].data
    CATALOG = h[3].data

    #Rebin data
    skewer_rows = utils.merge_cells(skewer_rows,N_merge)
    IVAR_rows = (utils.merge_cells(IVAR_rows,N_merge)==1).astype('float32')
    LOGLAM_MAP = np.log10(utils.merge_cells(10**LOGLAM_MAP,N_merge))

    #Filter out QSOs with less than the min number of relevant cells.
    relevant_QSOs = (np.sum(IVAR_rows,axis=1)>min_number_cells)
    skewer_rows = skewer_rows[relevant_QSOs,:]
    IVAR_rows = IVAR_rows[relevant_QSOs,:]
    CATALOG = CATALOG[relevant_QSOs]

    #Exponentiate to flux.
    skewer_rows = np.exp(-skewer_rows)

    #Divite by mean flux.
    s = fits.open(statistics_file)
    mean_F = s[1].data['F_MEAN']
    mean_F_z = s[1].data['z']
    mean_F = utils.merge_cells(mean_F,N_merge)
    skewer_rows /= mean_F
    s.close()

    #Reconstruct the non-delta HDUs.
    hdu_deltas_new = fits.PrimaryHDU(data=skewer_rows.T,header=header)
    hdu_iv_new = fits.ImageHDU(data=IVAR_rows.T,header=hdu_iv.header,name='IV')
    hdu_LOGLAM_MAP_new = fits.ImageHDU(data=LOGLAM_MAP,header=hdu_LOGLAM_MAP.header,name='LOGLAM_MAP')
    hdu_CATALOG_new = fits.BinTableHDU(CATALOG,header=hdu_CATALOG.header,name='CATALOG')

    hdulist = fits.HDUList([hdu_deltas_new, hdu_iv_new, hdu_LOGLAM_MAP_new, hdu_CATALOG_new])
    out_filepath = utils.get_file_name(dirname,'picca-flux-rebin-{}-new'.format(N_merge),16,pixel)
    hdulist.writeto(out_filepath,overwrite=overwrite)
    hdulist.close()
    h.close()

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

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = N_processes)
    results = []
    start_time = time.time()

    for pixel in pixels:
        pool.apply_async(convert_picca_tau_to_flux,(pixel),callback=log_result,error_callback=log_error)

    pool.close()
    pool.join()
