#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool
import time
import process_functions as functions

lya = 1215.67

basedir = '/Users/jfarr/Projects/repixelise/test_output/test_multi'
#basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_revamp/process_output_4096_32'
basedir = '/Users/James/Projects/test_data/pixel_0'

N_side = 8
pixels = list(range(0,1))

N_bins = 1000
lambda_lower = 800.0
IVAR_cutoff = lya

sum_delta = np.zeros(N_bins)
sum_delta_squared = np.zeros(N_bins)
N_contributions = np.zeros(N_bins)

binned_mean_delta = np.zeros(N_bins)
binned_var_delta = np.zeros(N_bins)

bins = np.linspace(lambda_lower,IVAR_cutoff,N_bins+1)
binned_lambdas = np.zeros(N_bins)
for i in range(N_bins):
    binned_lambdas[i] = (bins[i]+bins[i+1])/2

#For each pixel, get the files and MOCKIDs (HDU 1, 'MOCKID')
#Make a cross=pixel list
#Divide it up
#Make tasks (pixel,THING_ID list)
#Change c_p_c to run per task

def compute_pixel_contribution(pixel):

    sum_delta_pixel = np.zeros(N_bins)
    sum_delta_squared_pixel = np.zeros(N_bins)
    N_contributions_pixel = np.zeros(N_bins)

    #print('looking at pixel {}'.format(pixel))

    pixel_100 = pixel//100

    colore_filename = '{}/{}/{}/gaussian-colore-{}-{}.fits'.format(basedir,pixel_100,pixel,N_side,pixel)

    h = fits.open(colore_filename)

    DELTA_rows = h[2].data
    Z_QSO = h[1].data['Z_COSMO']
    Z = h[4].data['Z']
    LOGLAM_MAP = np.log10(lya*(1+Z))
    N_qso = DELTA_rows.shape[0]

    h.close()

    IVAR_rows = functions.make_IVAR_rows(IVAR_cutoff,Z_QSO,LOGLAM_MAP)

    LAMBDA_R_rows = np.zeros(DELTA_rows.shape)
    binned_LAMBDA_R_rows = np.zeros(DELTA_rows.shape)

    #print('calculating the lambda_R bins')
    lya_lambdas = 10**LOGLAM_MAP
    for i in range(N_qso):
        LAMBDA_R_rows[i,:] = ((lya_lambdas)/(1+Z_QSO[i]))*(IVAR_rows[i,:])

    binned_LAMBDA_R_rows = np.digitize(LAMBDA_R_rows,bins) - 1

    for n in range(N_bins):

        #print('looking at bin {} ({:4.1f} to {:4.1f}A)'.format(n,bins[n],bins[n+1]),end='\r')

        bin_coordinates_initial = np.where(binned_LAMBDA_R_rows == n)
        bin_coordinates_i = bin_coordinates_initial[0]
        bin_coordinates_j = bin_coordinates_initial[1]
        bin_coordinates = zip(bin_coordinates_i,bin_coordinates_j)

        sum_delta_bin = 0
        sum_delta_squared_bin = 0
        N_contributions_bin = 0
        for coordinate_pair in bin_coordinates:
            i = coordinate_pair[0]
            j = coordinate_pair[1]
            sum_delta_pixel[n] += DELTA_rows[i,j]
            sum_delta_squared_pixel[n] += (DELTA_rows[i,j])**2
            N_contributions_pixel[n] += 1


    return [pixel,sum_delta_pixel,sum_delta_squared_pixel,N_contributions_pixel]

#Get the number of processes.
N_processes = int(sys.argv[1])

#Define a progress-tracking function.
def log_result(retval):
    pixel = retval[0]
    #sum_delta_pixel = retval[1]
    #sum_delta_squared_pixel = retval[2]
    #N_contributions_pixel = retval[3]

    #sum_delta += sum_delta_pixel
    #sum_delta_squared += sum_delta_squared_pixel
    #N_contributions += N_contributions_pixel

    results.append(retval)
    N_complete = len(results)
    N_pixels = len(pixels)
    N_pixels_digits = int(np.log10(N_pixels)) + 1

    N_chunks = 20
    N_chunks_complete = int((N_complete*N_chunks) // (N_pixels))
    block_char = '-'
    progress_bar = '|' + block_char*N_chunks_complete + ' '*(N_chunks-N_chunks_complete) + '|'

    current_time = time.time()
    time_elapsed = current_time - start_time
    estimated_time_remaining = (time_elapsed)*((N_pixels-N_complete)/N_complete)

    print(' -> current progress: {} {:4d} of {:4d} complete ({:3.0%}), {:4.0f}s elapsed, ~{:4.0f}s remaining'.format(progress_bar,N_complete,N_pixels,N_complete/N_pixels,time_elapsed,estimated_time_remaining),end="\r")

#Define an error-tracking function.
def log_error(retval):
    print('error',retval)

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = N_processes)
    results = []

    start_time = time.time()

    print('calculating...')
    for pixel in pixels:
        pool.apply_async(compute_pixel_contribution,[pixel],callback=log_result,error_callback=log_error)

    pool.close()
    pool.join()

for result in results:
    sum_delta += result[1]
    sum_delta_squared += result[2]
    N_contributions += result[3]

for n in range(N_bins):
    if N_contributions[n] > 0:
        binned_mean_delta[n] = sum_delta[n]/N_contributions[n]
        binned_var_delta[n] = sum_delta_squared[n]/N_contributions[n] - (binned_mean_delta[n])**2

print(' ')


plt.figure()
#plt.errorbar(R_binned,xi,yerr=err_1,fmt='o')
plt.plot(binned_lambdas,binned_mean_delta)
plt.plot(binned_lambdas,binned_var_delta)
plt.plot(binned_lambdas,binned_mean_delta*(((binned_lambdas-IVAR_cutoff)/(np.ones(N_bins)*IVAR_cutoff-lambda_lower))**2)/(binned_mean_delta[0]))
plt.grid(True, which='both')
plt.axhline(y=0, color='k')
plt.savefig('xcf1d_{}_{}_with_squared.pdf'.format(pixels[0],pixels[-1]))

plt.figure()
#plt.errorbar(R_binned,xi,yerr=err_1,fmt='o')
plt.plot(binned_lambdas,binned_mean_delta)
plt.plot(binned_lambdas,binned_var_delta)
plt.grid(True, which='both')
plt.axhline(y=0, color='k')
plt.savefig('xcf1d_{}_{}.pdf'.format(pixels[0],pixels[-1]))

plt.figure()
#plt.errorbar(R_binned,xi,yerr=err_1,fmt='o')
plt.plot(binned_lambdas,binned_mean_delta)
limit = 1.5*max(binned_mean_delta[:int(N_bins/2.0)])
plt.ylim(-limit,limit)
plt.grid(True, which='both')
plt.axhline(y=0, color='k')
plt.savefig('xcf1d_{}_{}_zoom.pdf'.format(pixels[0],pixels[-1]))

plt.show()

"""
for pixel in pixels:

    sum_delta_pixel = np.zeros(N_bins)
    sum_delta_squared_pixel = np.zeros(N_bins)
    N_contributions_pixel = np.zeros(N_bins)

    print('looking at pixel {}'.format(pixel))

    pixel_100 = pixel//100

    colore_filename = '{}/{}/{}/gaussian-colore-{}-{}.fits'.format(basedir,pixel_100,pixel,N_side,pixel)

    h = fits.open(colore_filename)

    DELTA_rows = h[2].data
    Z_QSO = h[1].data['Z_COSMO']
    Z = h[4].data['Z']
    LOGLAM_MAP = np.log10(lya*(1+Z))
    N_qso = DELTA_rows.shape[0]

    h.close()

    IVAR_rows = functions.make_IVAR_rows(IVAR_cutoff,Z_QSO,LOGLAM_MAP)

    bins = np.linspace(lambda_lower,IVAR_cutoff,N_bins+1)
    binned_lambdas = np.zeros(N_bins)
    for i in range(N_bins):
        binned_lambdas[i] = (bins[i]+bins[i+1])/2

    LAMBDA_R_rows = np.zeros(DELTA_rows.shape)
    binned_LAMBDA_R_rows = np.zeros(DELTA_rows.shape)

    print('calculating the lambda_R bins')
    lya_lambdas = 10**LOGLAM_MAP
    for i in range(N_qso):
        LAMBDA_R_rows[i,:] = ((lya_lambdas)/(1+Z_QSO[i]))*(IVAR_rows[i,:])

    binned_LAMBDA_R_rows = np.digitize(LAMBDA_R_rows,bins) - 1

    for n in range(N_bins):

        print('looking at bin {} ({:4.1f} to {:4.1f}A)'.format(n,bins[n],bins[n+1]),end='\r')

        bin_coordinates_initial = np.where(binned_LAMBDA_R_rows == n)
        bin_coordinates_i = bin_coordinates_initial[0]
        bin_coordinates_j = bin_coordinates_initial[1]
        bin_coordinates = zip(bin_coordinates_i,bin_coordinates_j)

        sum_delta_bin = 0
        sum_delta_squared_bin = 0
        N_contributions_bin = 0
        for coordinate_pair in bin_coordinates:
            i = coordinate_pair[0]
            j = coordinate_pair[1]
            sum_delta[n] += DELTA_rows[i,j]
            sum_delta_squared[n] += (DELTA_rows[i,j])**2
            N_contributions[n] += 1

    print(' ')
"""
