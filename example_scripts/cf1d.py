#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import sys
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool
import time

# TODO: Change the retrival of R once the master file format has been modified.
# TODO: Convert the parallelisation to per-pixel
# TODO: Save the output data as a file

#basedir = '/Users/jfarr/Projects/repixelise/test_output/test_multi'
basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_revamp/test/lya1100'
basedir = sys.argv[2]

correl_quantity = 'gaussian'
N_side = 8
lambda_min = 3550
lya = 1215.67

rmax = 160.0
rmin = 0.0
nr = 40

N_bin_splits = 100

pixels = list(range(10,11))

skewers = []
IVAR_skewers = []
max_N_cells_picca = 0

#Get the skewers.
for pixel in pixels:

    pixel_100 = pixel//100

    #Determine the desired filenames.
    picca_filename = '{}/{}/{}/picca-{}-{}-{}.fits'.format(basedir,pixel_100,pixel,correl_quantity,N_side,pixel)
    colore_filename = '{}/{}/{}/physical-colore-{}-{}.fits'.format(basedir,pixel_100,pixel,N_side,pixel)

    #Open up the file/s.
    h_picca = fits.open(picca_filename)
    h_colore = fits.open(colore_filename)

    skewers += [h_picca[0].data[:,j] for j in range(h_picca[0].data.shape[1])]
    IVAR_skewers += [h_picca[1].data[:,j] for j in range(h_picca[1].data.shape[1])]

    max_N_cells_picca = max(max_N_cells_picca,skewers[-1].shape[0])

    #Re-make the full loglam map from colore data.
    colore_lya_lambdas = lya*(1+h_colore[4].data['Z'])

    #Use this along with the value of lambda_min to determine the section of 'R' that we want.
    first_index = np.argmax(colore_lya_lambdas > 3550)

    R_picca = h_colore[4].data['R'][first_index:first_index+max_N_cells_picca+1]

[h_picca[0].data[:,j] for j in range(h_picca[0].data.shape[0])]

#skewers = skewers[0:100]
N_skewers = len(skewers)
print('there are {} skewers in total.'.format(N_skewers))

mean_skewer = np.zeros(max_N_cells_picca)
for i in range(max_N_cells_picca):
    total = 0
    N = 0
    for skewer in skewers:
        if len(skewer) >= i:
            total += skewer[i]
            N += 1
    mean_skewer[i] = total/N


#Calculate the separations of pixel pairs. If separation > rmax, set to -1.
print('making binned separations')
bins = np.linspace(rmin,rmax,nr+1)
separations = np.zeros((max_N_cells_picca,max_N_cells_picca))
for i in range(max_N_cells_picca):
    for j in range(max_N_cells_picca):
        separation = abs(R_picca[i]-R_picca[j])
        if separation < rmax and separation >= rmin:
            separations[i,j] = separation
        elif separation >= rmax or separation < rmin:
            separations[i,j] = -1.0

#Put the seaparations into bins.
binned_separations = np.digitize(separations,bins)

#Subtract 1 from the binned_separations to start the bin numbering from 0.
#Also subtract an identity matrix as we do not want to include a pixel's correlation with itself.
binned_separations = binned_separations - 1 - np.identity(max_N_cells_picca)
R_binned = np.zeros(nr)
for i in range(nr):
    R_binned[i] = bins[i]+(bins[i+1]-bins[i])/2

#For each bin, make a list of the pixel pairs - coordinates - associated with that bin. Store a list of these lists.
bin_coordinates = []
for n in range(nr):
    coord_arrays = np.where(binned_separations==n)
    row_coords = coord_arrays[0]
    col_coords = coord_arrays[1]
    bin_coordinates += [list(zip(row_coords,col_coords))]

#Get the number of processes.
N_processes = int(sys.argv[1])

#Divide up the bins if desired. A larger number of jobs is helpful for judging progress as the program is running.
#N_bin_splits = N_processes//nr + 1
if N_bin_splits > 1:
    split_size = N_skewers//N_bin_splits
    splits = []
    for split_n in range(N_bin_splits):
        splits += [[split_size*split_n,split_size*(split_n+1)-1]]
    splits[N_bin_splits-1][1] = N_skewers-1
    tasks = []
    for bin_n in range(nr):
        for split in splits:
            tasks += [(bin_n,bin_coordinates[bin_n],skewers[split[0]:split[1]],IVAR_skewers[split[0]:split[1]])]
else:
    split_size = N_skewers
    tasks = [(bin_n,bin_coordinates[bin_n],skewers,IVAR_skewers) for bin_n in range(nr)]

print('divided job up into {} tasks, with ~{} skewers in each one.'.format(len(tasks),split_size))

#Define the worker function.
def get_bin_xi(bin_n,bin_coordinates,skewers,IVAR_skewers):
    del_squared_bin = 0
    N_contributions_bin = 0
    del_squared_dev_bin = 0
    N_skewers = len(skewers)

    #For each skewer, find the coordinates in the list which correspond to valid cells.
    #i.e. those pairs in which both cells have IVAR=1
    for skewer_n,skewer in enumerate(skewers):
        N_cells_skewer = skewer.shape[0]
        for coordinates in bin_coordinates:
            i = coordinates[0]
            j = coordinates[1]
            if i < N_cells_skewer and j < N_cells_skewer:
                if IVAR_skewers[skewer_n][i]*IVAR_skewers[skewer_n][j] != 0:

                    del_squared_pair = skewer[i]*skewer[j]
                    mean_pair = mean_skewer[i]*mean_skewer[j]

                    del_squared_bin += del_squared_pair
                    N_contributions_bin += 1
                    del_squared_dev_bin += del_squared_pair - mean_pair

    return [bin_n,N_contributions_bin,del_squared_bin,del_squared_dev_bin]

#Define a progress-tracking function.
def log_result(retval):
    bin_n = retval[0]
    N_contributions_bin = retval[1]
    del_squared_bin = retval[2]
    del_squared_dev_bin = retval[3]

    N_contributions[bin_n] += N_contributions_bin
    del_squared[bin_n] += del_squared_bin
    del_squared_dev[bin_n] += del_squared_dev_bin

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

    print(' -> current progress: {} {:4d} of {:4d} complete ({:3.0%}), {:4.0f}s elapsed, ~{:4.0f}s remaining'.format(progress_bar,N_complete,N_tasks,N_complete/N_tasks,time_elapsed,estimated_time_remaining),end="\r")

#Define an error-tracking function.
def log_error(retval):
    print('error',retval)

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = N_processes)
    results = []
    del_squared = np.zeros(nr)
    N_contributions = np.zeros(nr)
    del_squared_dev = np.zeros(nr)
    start_time = time.time()

    print('calculating xi')
    for task in tasks:
        pool.apply_async(get_bin_xi,task,callback=log_result,error_callback=log_error)

    #xi = pool.starmap(get_bin_xi,tasks)
    #print('making xi took {}s'.format(time.time()-start_time))

    pool.close()
    pool.join()

print(' ')

xi = del_squared/N_contributions
cov = del_squared_dev/N_contributions
err = np.sqrt(cov)

plt.figure()
plt.plot(R_binned,xi)
#plt.plot(R_binned,(mean*np.ones(nr)))
plt.xlabel('r [Mpc/h]')
plt.ylabel('xi(r)')
plt.grid(True, which='both')
plt.savefig('xi_{}_{}.pdf'.format(pixels[0],pixels[-1]))

plt.figure()
plt.plot(R_binned,xi*(R_binned**2))
#plt.plot(R_binned,(mean*R_binned**2))
plt.xlabel('r [Mpc/h]')
plt.ylabel('r^2 xi(r)')
plt.grid(True, which='both')
plt.savefig('xir2_{}_{}.pdf'.format(pixels[0],pixels[-1]))

plt.figure()
plt.errorbar(R_binned,xi*(R_binned**2),yerr=err*(R_binned**2),fmt='o')
#plt.plot(R_binned,(mean*R_binned**2))
plt.xlabel('r [Mpc/h]')
plt.ylabel('r^2 xi(r)')
plt.grid(True, which='both')
plt.savefig('xir2_{}_{}.pdf'.format(pixels[0],pixels[-1]))

#plt.show()
