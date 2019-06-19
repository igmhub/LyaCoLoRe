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

from lyacolore import Pk1D, tuning, utils

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

parser.add_argument('--z-width', type = float, default = 0.2, required=False,
                    help = 'width of z bins to use')

parser.add_argument('--file-type', type = str, default = 'flux', required=False,
                    help = 'type of file to measure from')

parser.add_argument('--units', type = str, default = 'km/s', required=False,
                    help = 'choose \"km/s\" or \"Mpc/h\"')

parser.add_argument('--show-plot', action="store_true", default = False, required=False,
                    help = 'do we want to show the plot or just save')

parser.add_argument('--save-data', action="store_true", default = False, required=False,
                    help = 'do we want to save the data')

parser.add_argument('--smoothing-radius', type = float, default = 25.0, required=False,
                    help = 'gaussian smoothing radius to account for')

parser.add_argument('--overwrite', action="store_true", default = False, required=False,
                    help = 'overwrite existing files')

parser.add_argument('--compressed-input', action="store_true", default = False, required=False,
                    help = 'compress output files to .fits.gz')

parser.add_argument('--k-min-plot', type = float, default = 0.001, required=False,
                    help = 'min value of k to plot')

parser.add_argument('--k-max-plot', type = float, default = 0.02, required=False,
                    help = 'max value of k to plot')

parser.add_argument('--N-k-values', type = int, default = None, required=False,
                    help = 'number of values of k to plot')

parser.add_argument('--Pk1D-data', type = str, default = None, required=False,
                    help = 'previously calculated Pk1d')

parser.add_argument('--k-plot-power', type = float, default = 0.0, required=False,
                    help = 'power of k by which to multiply Pk1D')

parser.add_argument('--add-error-areas', action="store_true", default = False, required=False,
                    help = 'plot errors in the form of a filled area')

parser.add_argument('--add-error-bars', action="store_true", default = False, required=False,
                    help = 'plot errors in the form of error bars')

parser.add_argument('--add-model-interval-areas', type = float, default = None, required=False,
                    help = 'plot areas of +/- %age on model values')

parser.add_argument('--add-model-interval-lines', type = float, default = None, required=False,
                    help = 'plot lines of +/- %age on model values')

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
z_values = np.array(args.z_values)
z_width = args.z_width
file_type = args.file_type
units = args.units
show_plot = args.show_plot
save_data = args.save_data
smoothing_radius = args.smoothing_radius
overwrite = args.overwrite

# TODO: print to confirm the arguments. e.g. "DLAs will be added"

if np.log2(N_side)-int(np.log2(N_side)) != 0:
    print('nside must be a power of 2!')
else:
    N_pix = 12*N_side**2

colours = ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'] * 2
fontsize = 16
plotsize = (12, 5)
dpi = 80

#Set style options everywhere.
#plt.rc('text', usetex=True)
plt.rc('font', size=fontsize)

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
Get P1D
"""

Pk1D_results = {}
if args.Pk1D_data is None:
    print('Getting Pk1D from each pixel...')

    #Get z and R along the skewers.
    m = fits.open(base_dir+'/master.fits')
    if 'colorecell' in file_type:
        z = m['COSMO_COL'].data['Z']
        R = m['COSMO_COL'].data['R']
    else:
        z = m['COSMO_EXP'].data['Z']
        R = m['COSMO_EXP'].data['R']
    m.close()

    #Determine if we're looking at the Gaussian skewers.
    gaussian = ('colorecell' in file_type)

    dr_hMpc = (R[-1] - R[0])/(R.shape[0] - 1)

    #Function to get deltas and ivar from each pixel.
    def get_pixel_P1D(pixel):
        dirname = utils.get_dir_name(base_dir,pixel)
        filename = utils.get_file_name(dirname,'picca-'+file_type,N_side,pixel,compressed=args.compressed_input)
        h = fits.open(filename)
        delta_rows = h[0].data.T
        ivar_rows = h[1].data.T
        z = 10**(h[2].data)/lya - 1
        Pk1D_results = {}
        for z_value in z_values:
            k, Pk, var = Pk1D.get_Pk1D(delta_rows,ivar_rows,dr_hMpc,z,z_value=z_value
                                    ,z_width=z_width,units=units,R1=smoothing_radius
                                    ,gaussian=gaussian)
            Pk1D_results[z_value] = {'k':k, 'Pk':Pk, 'var':var}
        return (pixel,Pk1D_results)

    tasks = [(pixel,) for pixel in pixels]

    #Run the multiprocessing pool
    if __name__ == '__main__':
        pool = Pool(processes = N_processes)
        results = []
        start_time = time.time()

        for task in tasks:
            pool.apply_async(get_pixel_P1D,task,callback=log_result,error_callback=log_error)

        pool.close()
        pool.join()

    #Combine the results from each pixel.
    for z_value in z_values:
        #This assumes each pixel has the same weight. Ideally would have different weights for each k-mode from each pixel.
        N = 1
        P_sum = results[0][1][z_value]['Pk']
        P2_sum = results[0][1][z_value]['Pk']**2
        for result in results[1:]:
            N += 1
            P_sum += result[1][z_value]['Pk']
            P2_sum += result[1][z_value]['Pk']**2
        Pk1D_zval_results = {'k': result[1][z_value]['k'],
                            'Pk': P_sum/N,
                            'var': (P2_sum/N - (P_sum/N)**2)}
        Pk1D_results[z_value] = Pk1D_zval_results
else:
    print('Loading pre-computed 1D power spectrum...')
    Pk1D_data = fits.open(args.Pk1D_data)
    N_pixels = Pk1D_data[0].header['N_pixels']
    for hdu in Pk1D_data[1:]:
        z_value = float(hdu.name)
        if z_value in z_values:
            Pk1D_results[float(hdu.name)] = {'k':hdu.data['k'], 'Pk':hdu.data['Pk'], 'var':hdu.data['var']}

################################################################################
"""
Save the data.
"""

#Save the data.
def save_P1D_values(Pk1D_results):

    header = fits.Header()
    header['units'] = units
    header['N_side'] = N_side
    header['N_pixels'] = N_pixels

    prihdu = fits.PrimaryHDU(header=header)
    hdus = [prihdu]

    for key in Pk1D_results.keys():

        #Extract the data from the results dictionary.
        k = Pk1D_results[key]['k']
        Pk = Pk1D_results[key]['Pk']
        var = Pk1D_results[key]['var']

        dtype = [('k', 'f4'), ('Pk', 'f4'), ('var', 'f4')]
        data = np.array(list(zip(k,Pk,var)),dtype=dtype)
        hdu = fits.BinTableHDU.from_columns(data,header=header,name=str(key))
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
Make a plot.
"""

#Plot the power spectra
def plot_P1D_values(Pk1D_results,show_plot=True):
    fig = plt.figure(figsize=(8, 5), dpi= 80, facecolor='w', edgecolor='k')

    #Extract the z values and sort them.
    z_values = np.array([key for key in Pk1D_results.keys()])
    z_values.sort()

    #Set up to record the max/min values of P so we can set axis scales well.
    plot_model_min = 10**10.
    plot_model_max = 0.

    #For each z value, plot Pk1D and a model line.
    for i,z_value in enumerate(z_values):

        #Extract the data from the results dictionary.
        k = Pk1D_results[z_value]['k']
        Pk = Pk1D_results[z_value]['Pk']
        var = Pk1D_results[z_value]['var']
        colour = colours[i]

        #Get the relevant k values for this plot.
        k_rel = (k>args.k_min_plot) * (k<args.k_max_plot)

        #Plot the result and a "model" comparison +/-10%
        #plt.loglog(k,Pk,label='z={}'.format(z_value),c=colour)
        to_plot_data = Pk * (k ** args.k_plot_power)
        to_plot_data_err = np.sqrt(var/N_pixels) * (k ** args.k_plot_power)
        if args.add_error_bars:
            if args.N_k_values is not None:
                #Extract the k values we want to plot.
                N_rel = np.sum(k_rel)
                N_low = np.sum(k<args.k_min_plot)
                ind_orig = np.arange(N_low,N_low+N_rel,1)
                spacing = N_rel // args.N_k_values
                ind = ind_orig[(ind_orig % spacing) == 0]
                print(k.shape,Pk.shape,ind.shape)
                plt.errorbar(k[ind],to_plot_data[ind],yerr=to_plot_data_err[ind],label='z={}'.format(z_value),c=colour,fmt='o')
            else:
                plt.errorbar(k,to_plot_data,yerr=to_plot_data_err,label='z={}'.format(z_value),c=colour,fmt='o')
        elif args.add_error_areas:
            plt.plot(k,to_plot_data,label='z={}'.format(z_value),c=colour)
            plt.fill_between(k,to_plot_data+to_plot_data_err,to_plot_data-to_plot_data_err,color=colour,alpha=0.2)
        else:
            plt.plot(k,to_plot_data,label='z={}'.format(z_value),c=colour)
        plt.semilogx()
        plt.semilogy()

        if 'flux' in file_type:
            model_Pk_kms = tuning.P1D_z_kms_PD2013(z_value,k)
            to_plot_model = model_Pk_kms * (k ** args.k_plot_power)
            plt.plot(k,to_plot_model,c=colour,linestyle='--')
            if args.add_model_interval_areas is not None:
                to_plot_lower = to_plot_model * (1 - args.add_model_interval_areas)
                to_plot_upper = to_plot_model * (1 + args.add_model_interval_areas)
                plt.fill_between(k,to_plot_upper,to_plot_lower,color=colour,alpha=0.2)
            elif args.add_model_interval_lines is not None:
                to_plot_lower = to_plot_model * (1 - args.add_model_interval_lines)
                to_plot_upper = to_plot_model * (1 + args.add_model_interval_lines)
                plt.plot(k,to_plot_upper,color=colour,linestyle=':')
                plt.plot(k,to_plot_lower,color=colour,linestyle=':')
            plot_model_min = np.minimum(plot_model_min,np.min(to_plot_model[k_rel]))
            plot_model_max = np.maximum(plot_model_max,np.max(to_plot_model[k_rel]))

        #plt.errorbar(k,Pk,yerr=np.sqrt(var),label='z={}'.format(key),marker='o',c=colour)
        #model_Pk_kms = tuning.P1D_z_kms_PD2013(key,k)
        #plt.plot(k,model_Pk_kms,label='DR9 z={}'.format(key),c=colour,linestyle=':')
        #plt.plot(k,model_Pk_kms*0.9,color=[0.5,0.5,0.5],alpha=0.5)
        #plt.plot(k,model_Pk_kms*1.1,color=[0.5,0.5,0.5],alpha=0.5)

    #plt.semilogy()
    #plt.semilogx()

    space = 0.1
    ylim_lower = plot_model_min * (1 - space)
    ylim_upper = plot_model_max * (1 + space)
    plt.xlim(args.k_min_plot,args.k_max_plot)
    plt.ylim(ylim_lower,ylim_upper)

    plt.legend()
    #plt.grid()
    #plt.subplots_adjust(right=0.65)
    #ax.legend(loc='upper left',bbox_to_anchor= (1.01, 1.0))
    ylabel = r'$P_{1D}$'
    if args.k_plot_power == 1 :
        ylabel = r'$k P_{{1D}}(k)\ /\ (kms^{{-1}})^{{-{:d}}}$'.format(int(args.k_plot_power))
    elif args.k_plot_power > 1 :
        ylabel = r'$k^{:d} P_{{1D}}(k)\ /\ (kms^{{-1}})^{{-{:d}}}$'.format(int(args.k_plot_power),int(args.k_plot_power))
    plt.ylabel(ylabel)
    if units == 'km/s':
        plt.xlabel(r'$k\ /\ (kms^{-1})^{-1}$')
    elif units == 'Mpc/h':
        plt.xlabel(r'$k\ /\ (Mpch^{-1})^{-1}$')
    filename = 'Pk1D_{}_{}.pdf'.format(file_type,N_pixels)
    plt.savefig('Pk1D_k{}_vs_k_{}.pdf'.format(int(args.k_plot_power),file_type))
    if show_plot:
        plt.show()
    return

plot_P1D_values(Pk1D_results,show_plot=show_plot)

################################################################################
"""
Celebrate!
"""
