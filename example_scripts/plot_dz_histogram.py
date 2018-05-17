import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
from multiprocessing import Pool
import multiprocessing
import sys
import time

lya = 1215.67
N_processes = int(sys.argv[1])

basedir = 'example_data/update_160518/'
basedir = '/Users/James/Projects/test_data/process_output_G_hZsmooth_4096_32_sr2.0_bm1_biasG18_picos_nside16/'
basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/test/'

files = glob.glob(basedir + '/*/*/gaussian-colore-16*.fits')

N_files = 1000

if N_files < len(files):
    files = files[:N_files]

N_bins = 100
z_min = 2.0
z_max = 2.2

rest_frame_cutoff = 1150.0 #Å

dz_upper = 0.01
dz_lower = -0.01
dz_range = (dz_lower,dz_upper)

QSO_hist = np.zeros(N_bins)
vel_skw_hist = np.zeros(N_bins)

N_QSO = 0
N_cells = 0

################################################################################

"""
Define the multiprocessing tracking functions
"""

#Define a progress-tracking function.
def log_result(retval):
    results.append(retval)
    N_complete = len(results)
    N_tasks = len(tasks)
    N_tasks_digits = int(np.log10(N_tasks)) + 1

    N_chunks = 20
    N_chunks_complete = int((N_complete*N_chunks) // (N_tasks))
    block_char = '-'
    progress_bar = '|' + block_char*N_chunks_complete + ' '*(N_chunks-N_chunks_complete) + '|'

    current_time = time.time()
    time_elapsed = current_time - start_time
    estimated_time_remaining = (time_elapsed)*((N_tasks-N_complete)/N_complete)
    print(' -> current progress: {} {:4d} of {:4d} complete ({:3.0%}), {:4.0f}s elapsed, ~{:5.0f}s remaining'.format(progress_bar,N_complete,N_tasks,N_complete/N_tasks,time_elapsed,estimated_time_remaining),end="\r")

    if len(results) == len(tasks):
        print('\nProcess complete!')

#Define an error-tracking function.
def log_error(retval):
    print('Error:',retval)

################################################################################

def make_file_histogram(file):

    h = fits.open(file)

    relevant_QSOs = [i for i in range(h[1].data.shape[0]) if h[1].data['Z_COSMO'][i] > z_min and h[1].data['Z_COSMO'][i] < z_max]

    file_QSO_hist,bins = np.histogram(h[1].data['DZ_RSD'][relevant_QSOs],bins=N_bins,range=dz_range)
    #N_QSO += len(relevant_QSOs)
    #QSO_hist += file_QSO_hist
    #print('done with QSOs')

    file_vel_skw_hist = np.zeros(N_bins)

    for i in range(h[1].data.shape[0]):
        #print('{}/{} skewers'.format(i,len(relevant_QSOs)),end='\r')

        Z_QSO = h[1].data['Z_COSMO'][i]
        Z = h[4].data['Z']

        last_cell = np.searchsorted(lya*(1+Z),rest_frame_cutoff*(1+Z_QSO))

        relevant_cells = [i for i in range(last_cell) if Z[i] > z_min and Z[i] < min(Z_QSO,z_max)]

        skewer_file_vel_skw_hist,bins = np.histogram(h[3].data[i,relevant_cells],bins=N_bins,range=dz_range)
        #N_cells += len(relevant_cells)
        file_vel_skw_hist += skewer_file_vel_skw_hist
    #print('done with vel skewers')

    h.close()

    return [bins,file_QSO_hist,file_vel_skw_hist]

tasks = files

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = N_processes)
    results = []
    start_time = time.time()

    for task in tasks:
        pool.apply_async(make_file_histogram,[task],callback=log_result,error_callback=log_error)

    pool.close()
    pool.join()

for result in results:
    bins = result[0]
    QSO_hist += result[1]
    vel_skw_hist += result[2]
    N_QSO += sum(result[1])
    N_cells += sum(result[2])

N_QSO = N_QSO.astype('int')
N_cells = N_cells.astype('int')
"""
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

"""

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
plt.savefig('dz_histogram_{}_{}.pdf'.format(N_files,rest_frame_cutoff))
plt.show()
