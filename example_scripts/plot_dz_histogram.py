import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob

basedir = 'example_data/update_160518/'
basedir = '/Users/jfarr/Projects/test_data/process_output_G_hZsmooth_4096_32_sr2.0_bm1_biasG18_picos_nside16/'
basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/test/'

files = glob.glob(basedir + '/*/*/gaussian-colore-16*.fits')

N_files = 100
N_QSO = 10000
N_cells = 10000

if N_files < len(files):
    files = files[:N_files]

N_bins = 100
z_min = 2.0
z_max = 2.2

QSO_DZs = []
vel_skw_values = []

for file in files:
    print(file)
    h = fits.open(file)

    if len(QSO_DZs)<N_QSO:
        relevant_QSOs = [i for i in range(h[1].data.shape[0]) if h[1].data['Z_COSMO'][i] > z_min and h[1].data['Z_COSMO'][i] < z_max]
        QSO_DZs += list(h[1].data['DZ_RSD'][relevant_QSOs])


    for relevant_QSO in relevant_QSOs:
        if len(vel_skw_values)<N_cells:
            relevant_cells = [i for i in range(h[4].data.shape[0]) if h[4].data['Z'][i] > z_min and h[4].data['Z'][i] < min(h[1].data['Z_COSMO'][relevant_QSO],z_max)]
            vel_skw_values += list(h[3].data[relevant_QSO,relevant_cells])

    h.close()

print(len(QSO_DZs),len(vel_skw_values))

plt.hist(QSO_DZs,bins=N_bins,normed=True,label='QSO dz values')
plt.hist(vel_skw_values,bins=N_bins,normed=True,label='vel skw values')
plt.legend()
plt.grid()
plt.show()
