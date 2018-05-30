import numpy as np
from astropy.io import fits

N_side = 16
N_files = 12*N_side**2
file_prefix = 'picca-gaussian-RSD'

"""
Formatting of statistics file:

statistics_dtype = [('N', 'f4')
    , ('GAUSSIAN_DELTA_MEAN', 'f4'), ('GAUSSIAN_DELTA_VAR', 'f4')
    , ('DENSITY_DELTA_MEAN', 'f4'), ('DENSITY_DELTA_VAR', 'f4')
    , ('F_MEAN', 'f4'), ('F_VAR', 'f4')
    , ('F_DELTA_MEAN', 'f4'), ('F_DELTA_VAR', 'f4')]
"""

statistics = fits.open('statistics.fits')
mean = statistics['GAUSSIAN_DELTA_MEAN']
N_cells = mean.shape[0]

for N in range(N_files):
    filepath = '{}/{}/'.format(N//100,N) + file_prefix + '-{}-{}.fits'.format(N_side,N)
    h = fits.open(filepath)
    GAUSSIAN_DELTA_rows = h[0].data.T

    renormalised_GAUSSIAN_DELTA_rows = np.zeros(GAUSSIAN_DELTA_rows.shape)
    for j in range(N_cells):
        if mean[j] != 0.0:
            renormalised_GAUSSIAN_DELTA_rows[:,j] = GAUSSIAN_DELTA_rows[:,j]/mean[j]

    header = h[0].header
    picca_0 = renormalised_GAUSSIAN_DELTA_rows.T
    hdu_DELTA = fits.PrimaryHDU(data=picca_0,header=header)

    hdu_iv = h[1]
    hdu_LOGLAM_MAP = h[2]
    hdu_CATALOG = h[3]
    h.close()

    #Combine the HDUs into and HDUlist and save as a new file. Close the HDUlist.
    hdulist = fits.HDUList([hdu_DELTA, hdu_iv, hdu_LOGLAM_MAP, hdu_CATALOG])

    new_filepath = '{}/{}/'.format(N//100,N) + file_prefix + '-renorm-{}-{}.fits'.format(N_side,N)
    hdulist.writeto(location+filename)
    hdulist.close()
