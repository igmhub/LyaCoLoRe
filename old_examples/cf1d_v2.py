#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool
import time

# TODO: Change the retrival of R once the master file format has been modified.
# TODO: Change the skewer getting system to do it on the fly

#basedir = '/Users/jfarr/Projects/repixelise/test_output/test_multi'
#basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_revamp/test/lya1100/'
basedir = sys.argv[2]

save_location = '/global/homes/j/jfarr/Projects/LyaCoLoRe/'

correl_quantity = 'gaussian'
N_side = 8
lambda_min = 3550
lya = 1215.67

rmax = 160.0
rmin = 0.0
nr = 40

N_bin_splits = 100

pixels = list(range(int(sys.argv[3]),int(sys.argv[4])))

skewers_list = []
IVAR_skewers_list = []

#Get the skewers
for pixel in pixels:
    pixel_100 = pixel//100
    picca_filename = '{}/{}/{}/picca-{}-{}-{}.fits'.format(basedir,pixel_100,pixel,correl_quantity,N_side,pixel)

    h_picca = fits.open(picca_filename)

    skewers_list += [h_picca[0].data[:,j] for j in range(h_picca[0].data.shape[1])]
    IVAR_skewers_list += [h_picca[1].data[:,j] for j in range(h_picca[1].data.shape[1])]

    h_picca.close()

skewers = np.asarray(skewers_list).T
IVAR_skewers = np.asarray(IVAR_skewers_list).T
N_cells = skewers.shape[0]

def get_skewers(pixel):
    pixel_100 = pixel//100
    picca_filename = '{}/{}/{}/picca-{}-{}-{}.fits'.format(basedir,pixel_100,pixel,correl_quantity,N_side,pixel)

    h_picca = fits.open(picca_filename)

    skewers = h_picca[0].data
    IVAR_skewers = h_picca[1].data

    h_picca.close()

    return skewers, IVAR_skewers

#Get the R vector.
# TODO: replace this with capability to retrieve R from a new master HDU extension
colore_filename = '{}/{}/{}/physical-colore-{}-{}.fits'.format(basedir,pixel_100,pixels[0],N_side,pixels[0])
h_colore = fits.open(colore_filename)

colore_lya_lambdas = lya*(1+h_colore[4].data['Z'])

first_index = np.argmax(colore_lya_lambdas > 3550)
R_picca = h_colore[4].data['R'][first_index:]

h_colore.close()

#skewers = skewers[:,0:100]
N_skewers = skewers.shape[1]
print('there are {} skewers in total, each with {} cells.'.format(N_skewers,N_cells))

#Calculate the separations of pixel pairs. If separation > rmax, set to -1.
print('making binned separations')
bins = np.linspace(rmin,rmax,nr+1)
separations = np.zeros((N_cells,N_cells))
for i in range(N_cells):
    for j in range(N_cells):
        separation = abs(R_picca[i]-R_picca[j])
        if separation < rmax and separation >= rmin:
            separations[i,j] = separation
        elif separation >= rmax or separation < rmin:
            separations[i,j] = -1.0

#Put the seaparations into bins.
binned_separations = np.digitize(separations,bins)

#Subtract 1 from the binned_separations to start the bin numbering from 0.
#Also subtract an identity matrix as we do not want to include a pixel's correlation with itself.
binned_separations = binned_separations - 1 - np.identity(N_cells)
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
            tasks += [(bin_n,bin_coordinates[bin_n],skewers[:,split[0]:split[1]],IVAR_skewers[:,split[0]:split[1]])]
else:
    split_size = N_skewers
    tasks = [(bin_n,bin_coordinates[bin_n],skewers,IVAR_skewers) for bin_n in range(nr)]

print('divided job up into {} tasks, with ~{} skewers in each one.'.format(len(tasks),split_size))

#Define the worker function.
def get_xi(bin_n,bin_coordinates,skewers,IVAR_skewers):
    del_squared_chunk = 0
    N_contributions_chunk = 0
    N_skewers = skewers.shape[1]

    #For each skewer, find the coordinates in the list which correspond to valid cells.
    #i.e. those pairs in which both cells have IVAR=1
    for skewer_n in range(N_skewers):
        skewer = skewers[:,skewer_n]
        IVAR_skewer = IVAR_skewers[:,skewer_n]
        for coordinates in bin_coordinates:
            i = coordinates[0]
            j = coordinates[1]

            if IVAR_skewer[i]*IVAR_skewer[j] != 0:
                del_squared_chunk += skewer[i]*skewer[j]
                N_contributions_chunk += 1

    if N_contributions_chunk > 0:
        xi_chunk = del_squared_chunk/N_contributions_chunk
    else:
        xi_chunk = 0

    return [bin_n,N_contributions_chunk,xi_chunk,del_squared_chunk]

#Define a progress-tracking function.
def log_result(retval):
    bin_n = retval[0]
    N_contributions_bin = retval[1]
    xi_bin = retval[2]
    del_squared_bin = retval[3]

    N_contributions[bin_n] += N_contributions_bin
    del_squared[bin_n] += del_squared_bin

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
    N_contributions = np.zeros(nr)
    del_squared = np.zeros(nr)
    start_time = time.time()

    print('calculating xi')
    for task in tasks:
        pool.apply_async(get_xi,task,callback=log_result,error_callback=log_error)

    #xi = pool.starmap(get_bin_xi,tasks)
    #print('making xi took {}s'.format(time.time()-start_time))

    pool.close()
    pool.join()

print(' ')

xi_alt = del_squared/N_contributions

# TODO: parallelise this
#Making variance to get the error bar
print('estimating errors')
xi = np.zeros(nr)
mean_xi_squared = np.zeros(nr)
var_xi = np.zeros(nr)

for bin_n in range(nr):

    weight_sum = 0

    for result in results:
        if result[0] == bin_n:

            result_N_contributions = result[1]

            if result_N_contributions != 0:
                weight = result_N_contributions/N_contributions[bin_n]
                result_xi = result[2]
            else:
                weight = 0
                result_xi = 0

            weighted_result_xi = weight*result_xi
            xi[bin_n] += weighted_result_xi
            mean_xi_squared[bin_n] += weight*(result_xi**2)

            weight_sum += weight

    #print(var_xi[bin_n])
    var_xi[bin_n] = mean_xi_squared[bin_n] - (xi[bin_n])**2

err = np.sqrt(var_xi)

plt.figure()
plt.plot(R_binned,xi)
#plt.plot(R_binned,xi_alt)
#plt.plot(R_binned,(mean*np.ones(nr)))
plt.xlabel('r [Mpc/h]')
plt.ylabel('xi(r)')
plt.grid(True, which='both')
plt.savefig(save_location+'xi_v2_{}_{}_{}.pdf'.format(pixels[0],pixels[-1],basedir[-7:-2]))

plt.figure()
plt.errorbar(R_binned,xi,yerr=err,fmt='o')
#plt.plot(R_binned,(mean*R_binned**2))
plt.xlabel('r [Mpc/h]')
plt.ylabel('xi(r)')
plt.grid(True, which='both')
plt.savefig(save_location+'xi_v2_{}_{}_errors_{}.pdf'.format(pixels[0],pixels[-1],basedir[-7:-2]))

plt.figure()
plt.plot(R_binned,xi*(R_binned**2))
#plt.plot(R_binned,(mean*R_binned**2))
plt.xlabel('r [Mpc/h]')
plt.ylabel('r^2 xi(r)')
plt.grid(True, which='both')
plt.savefig(save_location+'xir2_v2_{}_{}_{}.pdf'.format(pixels[0],pixels[-1],basedir[-7:-2]))

plt.figure()
plt.errorbar(R_binned,xi*(R_binned**2),yerr=err*(R_binned**2),fmt='o')
#plt.plot(R_binned,(mean*R_binned**2))
plt.xlabel('r [Mpc/h]')
plt.ylabel('r^2 xi(r)')
plt.grid(True, which='both')
plt.savefig(save_location+'xir2_v2_{}_{}_errors_{}.pdf'.format(pixels[0],pixels[-1],basedir[-7:-2]))

file_data_list = []
for result in results:
    file_data_list += [(result[0],result[1],result[2],result[3])]

dtype = [('bin_n', '>f4'), ('N_contributions_chunk', '>f4'), ('xi_chunk', '>f4'), ('del_squared_chunk', '>f4')]
file_data = np.array(file_data_list,dtype=dtype)

prihdr = fits.Header()
prihdu = fits.PrimaryHDU(header=prihdr)
cols_xi = fits.ColDefs(file_data)
hdu_xi = fits.BinTableHDU.from_columns(cols_xi,name='XI DATA')

hdulist = fits.HDUList([prihdu, hdu_xi])
hdulist.writeto(save_location+'xi_data_v2_{}_{}_{}.fits'.format(pixels[0],pixels[-1],basedir[-7:-2]))
hdulist.close

print('files saved to {}!'.format(save_location))
#plt.show()
