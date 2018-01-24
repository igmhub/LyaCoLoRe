import numpy as np
from astropy.io import fits
import sys
import matplotlib.pyplot as plt
from multiprocessing import Pool
import time


basedir = '/Users/jfarr/Projects/repixelise/test_output/test_multi'
basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_revamp/process_output_4096_32'
correl_quantity = 'gaussian'
N_side = 8
lambda_min = 3550
lya = 1215.67

rmax = 160.0
rmin = 0.0
nr = 40

pixels = list(range(10,20))

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

#skewers = skewers[-1000:]
N_skewers = len(skewers)
print('there are {} skewers in total.'.format(N_skewers))


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

"""
OBSOLETE

#Produce a list of nr masks, one for each r bin.
reduced_separations = binned_separations + 1 + np.identity(max_N_cells_picca)
bin_masks = []
for n in range(nr,0,-1):
    mask = reduced_separations//n
    bin_masks += [mask]
    reduced_separations = reduced_separations - n*mask

bin_masks[-1] = bin_masks[-1] - np.identity(bin_masks[-1].shape[0])
bin_masks.reverse()
"""

#Convert each of these masks into a list of coordinates. Store a list of these lists.
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
N_bin_splits = 1
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

                    del_squared_bin += skewer[i]*skewer[j]
                    N_contributions += 1

    return [bin_n,N_contributions_bin,del_squared_bin]

#Define a progress-tracking function.
def log_result(retval):
    bin_n = retval[0]
    N_contributions_bin = retval[1]
    del_squared_bin = retval[2]

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
    del_squared = np.zeros(nr)
    N_cells = np.zeros(nr)
    start_time = time.time()

    print('calculating xi')
    for task in tasks:
        a=[]
        pool.apply_async(get_bin_xi,task,callback=log_result,error_callback=log_error)

    #xi = pool.starmap(get_bin_xi,tasks)
    #print('making xi took {}s'.format(time.time()-start_time))

    pool.close()
    pool.join()

print(' ')

xi = del_squared/N_cells

mean_delta = np.zeros(max_N_cells_picca)
for i in range(max_N_cells_picca):
    count = 0
    sum = 0
    for n in range(N_skewers):
        if len(skewers[n]) > i and IVAR_skewers[n][i] == 1:
            sum += skewers[n][i]
            count += 1
    if count != 0:
        mean_delta[i] = sum/count
    else:
        mean_delta[i] = 0

mean = np.mean(mean_delta)
print('mean delta is {}'.format(mean))

plt.figure()
plt.plot(R_binned,xi)
#plt.plot(R_binned,(mean*np.ones(nr)))

plt.figure()
plt.plot(R_binned,xi*(R_binned**2))
#plt.plot(R_binned,(mean*R_binned**2))

plt.show()
