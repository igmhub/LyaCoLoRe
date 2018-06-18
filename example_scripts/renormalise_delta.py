import numpy as np
from astropy.io import fits
import sys

N_side = 16
if len(sys.argv) >= 1:
    N_files = int(sys.argv[1])
else:
    N_files = 12*N_side**2

file_prefix = 'picca-gaussian-noRSD'
basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/test/'

"""
Formatting of statistics file:

statistics_dtype = [('N', 'f4')
    , ('GAUSSIAN_DELTA_MEAN', 'f4'), ('GAUSSIAN_DELTA_VAR', 'f4')
    , ('DENSITY_DELTA_MEAN', 'f4'), ('DENSITY_DELTA_VAR', 'f4')
    , ('F_MEAN', 'f4'), ('F_VAR', 'f4')
    , ('F_DELTA_MEAN', 'f4'), ('F_DELTA_VAR', 'f4')]
"""

"""
statistics = fits.open(basedir + 'statistics.fits')
mean = statistics[1].data['GAUSSIAN_DELTA_MEAN']
"""

N_skewers = 0
mean = 0

N_skewers = []
mean = []

for N in range(N_files):
    print('calculating mean from file {} of {}'.format(N+1,N_files),end='\r')

    filepath = basedir + '/{}/{}/'.format(N//100,N) + file_prefix + '-{}-{}.fits'.format(N_side,N)
    h = fits.open(filepath)

    data_rows = h[0].data.T
    weights_rows = h[1].data.T

    N_skewers_file = np.sum(weights_rows,axis=0)
    N_skewers += [N_skewers_file]
    mean += [np.sum(data_rows*weights_rows,axis=0)]

    h.close()

mean = np.sum(mean,axis=0)
N_skewers = np.sum(N_skewers,axis=0)

N_cells = mean.shape[0]

for j in range(N_cells):
    if N_skewers[j] != 0:
        mean[j] /= N_skewers[j]

print('\n')
print(mean)

for N in range(N_files):
    filepath = basedir + '/{}/{}/'.format(N//100,N) + file_prefix + '-{}-{}.fits'.format(N_side,N)
    h = fits.open(filepath)
    GAUSSIAN_DELTA_rows = h[0].data.T

    renormalised_GAUSSIAN_DELTA_rows = np.zeros(GAUSSIAN_DELTA_rows.shape)
    for j in range(N_cells):
        if mean[j] != 0.0:
            renormalised_GAUSSIAN_DELTA_rows[:,j] = (1 + GAUSSIAN_DELTA_rows[:,j])/(1 + mean[j]) - 1

    header = h[0].header
    picca_0 = renormalised_GAUSSIAN_DELTA_rows.T
    hdu_DELTA = fits.PrimaryHDU(data=picca_0,header=header)

    hdu_iv = h[1]
    hdu_LOGLAM_MAP = h[2]
    hdu_CATALOG = h[3]

    #Combine the HDUs into and HDUlist and save as a new file. Close the HDUlist.
    hdulist = fits.HDUList([hdu_DELTA, hdu_iv, hdu_LOGLAM_MAP, hdu_CATALOG])

    new_filepath = basedir + '/{}/{}/'.format(N//100,N) + file_prefix + '-renorm-{}-{}.fits'.format(N_side,N)
    hdulist.writeto(new_filepath,overwrite=True)
    hdulist.close()
    h.close()

    print('{}/{} complete'.format(N+1,N_files),end='\r')
