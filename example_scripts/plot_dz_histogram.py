import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob

basedir = 'example_data/update_160518/'
basedir = '/Users/jfarr/Projects/test_data/process_output_G_hZsmooth_4096_32_sr2.0_bm1_biasG18_picos_nside16/'
basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/test/'

files = glob.glob(basedir + '/*/*/gaussian-colore-16*.fits')

N_files = 250

if N_files < len(files):
    files = files[:N_files]

N_bins = 100
z_min = 2.0
z_max = 2.2

dz_upper = 0.01
dz_lower = -0.01
dz_range = (dz_lower,dz_upper)

QSO_hist = np.zeros(N_bins)
vel_skw_hist = np.zeros(N_bins)

N_QSO = 0
N_cells = 0

for k,file in enumerate(files):
    print('\nlooking at {} ({}/{})'.format(file,k+1,len(files)))
    h = fits.open(file)

    relevant_QSOs = [i for i in range(h[1].data.shape[0]) if h[1].data['Z_COSMO'][i] > z_min and h[1].data['Z_COSMO'][i] < z_max]
    file_QSO_hist,bins = np.histogram(h[1].data['DZ_RSD'][relevant_QSOs],bins=N_bins,range=dz_range)
    N_QSO += len(relevant_QSOs)
    QSO_hist += file_QSO_hist
    print('done with QSOs')

    for i,relevant_QSO in enumerate(relevant_QSOs):
        print('{}/{} skewers'.format(i,len(relevant_QSOs)),end='\r')
        relevant_cells = [i for i in range(h[4].data.shape[0]) if h[4].data['Z'][i] > z_min and h[4].data['Z'][i] < min(h[1].data['Z_COSMO'][relevant_QSO],z_max)]
        file_vel_skw_hist,bins = np.histogram(h[3].data[relevant_QSO,relevant_cells],bins=N_bins,range=dz_range)
        N_cells += len(relevant_cells)
        vel_skw_hist += file_vel_skw_hist
    print('done with vel skewers')

    h.close()

print('looked at {} QSO dz values'.format(N_QSO))
print('looked at {} vel skewer cells'.format(N_cells))

bin_centres = []
for i in range(len(bins)-1):
    bin_centres += [(bins[i]+bins[i+1])/2]

QSO_hist /= sum(QSO_hist)
vel_skw_hist /= sum(vel_skw_hist)

plt.plot(bin_centres,QSO_hist,label='QSO dz values')
plt.plot(bin_centres,vel_skw_hist,label='vel skw values')
plt.legend()
plt.grid()
plt.savefig('dz_histogram.pdf')
plt.show()
