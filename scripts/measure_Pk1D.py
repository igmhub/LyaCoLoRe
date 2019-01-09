#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from multiprocessing import Pool
import multiprocessing
import sys
import time
import os
import argparse

from pyacolore import Pk1D, tuning, utils

################################################################################

"""
Set up the file locations and filename structures.
Also define option preferences.
"""

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--base-dir', type = str, default = None, required=True,
                    help = 'base data directory')

parser.add_argument('--nproc', type = int, default = 1, required=False,
                    help = 'number of processes to use')

parser.add_argument('--nside', type = int, default = 16, required=False,
                    help = 'HEALPix nside for output files (must be 2^n)')

parser.add_argument('--pixels', type = int, default = None, required=False,
                    help = 'which pixel numbers to work on', nargs='*')

parser.add_argument('--z-values', type = float, default = [2.4], required=False,
                    help = 'which z values to measure at', nargs='*')

parser.add_argument('--z-width', type = int, default = 0.2, required=False,
                    help = 'width of z bins to use')

parser.add_argument('--file-type', type = str, default = 'flux', required=False,
                    help = 'type of file to measure from')

parser.add_argument('--units', type = str, default = 'km/s', required=False,
                    help = 'choose \"km/s\" or \"Mpc/h\"')

parser.add_argument('--show-plot', action="store_true", default = False, required=False,
                    help = 'do we want to show the plot or just save')

parser.add_argument('--save-data', action="store_true", default = False, required=False,
                    help = 'do we want to save the data')

################################################################################

print('setup arguments from parser')

args = parser.parse_args()

#Define global variables.
lya = utils.lya_rest

base_dir = args.base_dir
N_side = args.nside
pixels = args.pixels
if not pixels:
    pixels = list(range(12*N_side**2))
N_pixels = len(pixels)
N_processes = args.nproc
z_values = args.z_values
z_width = args.z_width
file_type = args.file_type
units = args.units
show_plot = args.show_plot
save_data = args.save_data

# TODO: print to confirm the arguments. e.g. "DLAs will be added"

if np.log2(N_side)-int(np.log2(N_side)) != 0:
    print('nside must be a power of 2!')
else:
    N_pix = 12*N_side**2

colours = ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'] * 2

################################################################################

"""
Define the multiprocessing tracking functions
"""

#Define a progress-tracking function.
def log_result(retval):

    results.append(retval)
    N_complete = len(results)
    N_tasks = len(tasks)

    utils.progress_bar(N_complete,N_tasks,start_time)

#Define an error-tracking function.
def log_error(retval):
    print('Error:',retval)

################################################################################
"""
Get the data
"""

print('Assembling the skewers data...')

#Get z and R along the skewers.
m = fits.open(base_dir+'/master.fits')
if 'colorecell' in file_type:
    z = m['COSMO_COL'].data['Z']
    R = m['COSMO_COL'].data['R']
else:
    z = m['COSMO_EXP'].data['Z']
    R = m['COSMO_EXP'].data['R']
m.close()

dr_hMpc = (R[-1] - R[0])/(R.shape[0] - 1)

#Function to get deltas and ivar from each pixel.
def get_pixel_data(pixel):
    dirname = utils.get_dir_name(base_dir,pixel)
    filename = utils.get_file_name(dirname,'picca-'+file_type,N_side,pixel)
    h = fits.open(filename)
    delta_rows = h[0].data.T
    ivar_rows = h[1].data.T
    return (delta_rows,ivar_rows)

tasks = [(pixel,) for pixel in pixels]

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = N_processes)
    results = []
    start_time = time.time()

    for task in tasks:
        pool.apply_async(get_pixel_data,task,callback=log_result,error_callback=log_error)

    pool.close()
    pool.join()

#Combine the results from each pixel.
delta_rows = results[0][0]
ivar_rows = results[0][1]
for result in results[1:]:
    delta_rows = np.append(delta_rows,result[0],axis=0)
    ivar_rows = np.append(ivar_rows,result[1],axis=0)

################################################################################
"""
Compute the P1D
"""

print('Computing the 1D power spectrum...')

Pk1D_results = {}
for i,z_value in enumerate(z_values):
    #print(z_value,delta_rows.shape,ivar_rows.shape,R.shape,z.shape)
    k, Pk, var = Pk1D.get_Pk1D(delta_rows,ivar_rows,dr_hMpc,z,z_value=z_value,z_width=z_width,units=units)
    Pk1D_results[z_value] = {'k':k, 'Pk':Pk, 'var':var}

################################################################################
"""
Make a plot.
"""

#Plot the power spectra
def plot_P1D_values(Pk1D_results,show_plot=True):
    plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')

    for key in Pk1D_results.keys():

        #Extract the data from the results dictionary.
        k = Pk1D_results[key]['k']
        Pk = Pk1D_results[key]['Pk']
        var = Pk1D_results[key]['var']
        colour = colours[np.searchsorted(z_values,key)]

        #Plot the result and a "model" comparison +/-10%
        plt.loglog(k,Pk,label='z={}'.format(key),c=colour)
        if 'flux' in file_type:
            model_Pk_kms = tuning.P1D_z_kms_PD2013(key,k)
            plt.loglog(k,model_Pk_kms,label='DR9 z={}'.format(key),c=colour,linestyle=':')
            plt.loglog(k,model_Pk_kms*0.9,color=[0.5,0.5,0.5],alpha=0.5)
            plt.loglog(k,model_Pk_kms*1.1,color=[0.5,0.5,0.5],alpha=0.5)

        #plt.errorbar(k,Pk,yerr=np.sqrt(var),label='z={}'.format(key),marker='o',c=colour)
        #model_Pk_kms = tuning.P1D_z_kms_PD2013(key,k)
        #plt.plot(k,model_Pk_kms,label='DR9 z={}'.format(key),c=colour,linestyle=':')
        #plt.plot(k,model_Pk_kms*0.9,color=[0.5,0.5,0.5],alpha=0.5)
        #plt.plot(k,model_Pk_kms*1.1,color=[0.5,0.5,0.5],alpha=0.5)

    #plt.semilogy()
    #plt.semilogx()

    #ylim_lower = min(model_Pk_kms) * 0.8
    #ylim_upper = max(model_Pk_kms) * 1.2
    #plt.ylim(ylim_lower,ylim_upper)

    plt.grid()
    plt.legend(fontsize=12)
    plt.ylabel(r'$P_{1D}$',fontsize=12)
    if units == 'km/s':
        plt.xlabel(r'$k\ /\ (kms^{-1})^{-1}$',fontsize=12)
    elif units == 'Mpc/h':
        plt.xlabel(r'$k\ /\ (Mpch^{-1})^{-1}$',fontsize=12)
    filename = 'Pk1D_{}_{}.pdf'.format(file_type,N_pixels)
    plt.savefig('Pk1D_{}.pdf'.format(file_type))
    if show_plot:
        plt.show()
    return

plot_P1D_values(Pk1D_results,show_plot=show_plot)

################################################################################
"""
Save the data.
"""

#Save the data.
def save_P1D_values(Pk1D_results):

    header = fits.Header()
    header['units'] = units
    header['N_side'] = N_side

    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)
    hdus = [prihdu]

    for key in Pk1D_results.keys():

        #Extract the data from the results dictionary.
        k = Pk1D_results[key]['k']
        Pk = Pk1D_results[key]['Pk']
        var = Pk1D_results[key]['var']

        dtype = [('k', 'f8'), ('Pk', 'f8'), ('var', 'f8')]
        data = np.array(list(zip(k,Pk,var)),dtype=dtype)
        hdu = fits.BinTableHDU.from_columns(data,header=header,name=key)
        hdus += [hdu]

    #Combine the HDUs into an HDUlist and save as a new file. Close the HDUlist.
    hdulist = fits.HDUList(hdus)
    filename = 'Pk1D_data_{}_{}.fits'.format(file_type,N_pixels)
    hdulist.writeto(filename,overwrite=overwrite)
    hdulist.close

    return

if save_data:
    save_P1D_values(Pk1D_results)

################################################################################
"""
Celebrate!
"""
