import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
from multiprocessing import Pool
import multiprocessing
import sys
import time

import utils

lya = utils.lya_rest

N_processes = int(sys.argv[1])

basedir = '../example_data/lya_skewers/'

files = glob.glob(basedir + '/*/*/gaussian-colore-16*.fits')

N_files = 3072

if N_files < len(files):
    files = files[:N_files]

N_bins = 100
z_min = 2.0
z_max = 2.2

rest_frame_cutoff = 1100. #Å

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

    #relevant_QSOs = [i for i in range(h[1].data.shape[0]) if h[1].data['Z_COSMO'][i] > z_min and h[1].data['Z_COSMO'][i] < z_max]
    #file_QSO_hist,bins = np.histogram(h[1].data['DZ_RSD'][relevant_QSOs],bins=N_bins,range=dz_range)
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

    return [bins,file_vel_skw_hist]

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


print('looking at master file to make QSO DZ_RSD histogram...')
master = fits.open(basedir + 'master.fits')
Z_QSO = master[1].data['Z_QSO_NO_RSD']
DZ_RSD = master[1].data['Z_QSO_RSD'] - master[1].data['Z_QSO_NO_RSD']

relevant_QSOs = [i for i in range(Z_QSO.shape[0]) if Z_QSO[i] > z_min and Z_QSO[i] < z_max]
QSO_hist,bins = np.histogram(DZ_RSD[relevant_QSOs],bins=N_bins,range=dz_range)

QSO_hist = np.array(QSO_hist,dtype='float')
N_QSO = len(relevant_QSOs)

for result in results:
    bins = result[0]
    #QSO_hist += result[1]
    vel_skw_hist += result[1]
    #N_QSO += sum(result[1])
    N_cells += sum(result[1])

#N_QSO = N_QSO.astype('int')
N_cells = N_cells.astype('int')

print('looked at {} QSO dz values'.format(N_QSO))
print('looked at {} vel skewer cells'.format(N_cells))

bin_centres = []
for i in range(len(bins)-1):
    bin_centres += [(bins[i]+bins[i+1])/2]

QSO_hist /= sum(QSO_hist)
vel_skw_hist /= sum(vel_skw_hist)

plt.plot(bin_centres,QSO_hist,label='QSO dz values')
plt.plot(bin_centres,vel_skw_hist,label='vel skw values')
plt.xlabel('dz')
plt.title('Histogram of RSD dz from QSO and skewer cells')
plt.legend()
plt.grid()
plt.savefig('dz_histogram_{}_{}.pdf'.format(N_files,rest_frame_cutoff))
plt.show()


