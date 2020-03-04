#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.ticker
from multiprocessing import Pool
import multiprocessing
import sys
import time
import os
import argparse
import matplotlib.patches as mpatches

from lyacolore import Pk1D, tuning, utils

################################################################################

"""
Set up the file locations and filename structures.
Also define option preferences.
"""

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--base-dir', type = str, default = None, required=False,
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

parser.add_argument('--quantity', type = str, default = 'flux', required=False,
                    help = 'quantity represented in file')

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

parser.add_argument('--add-vline', type = float, default = None, required=False,
                    help = 'add vertical line')

parser.add_argument('--fill-right-vline', action="store_true", default = False, required=False,
                    help = 'shade to the right of the vertical line')

parser.add_argument('--ysubs', type = float, default = None, required=False,
                    help = 'multiples of integer powers of 10 at which to place y ticks', nargs='*')

parser.add_argument('--save-suffix', type = str, default = None, required=False,
                    help = 'suffix to add onto save names for plots and data')

################################################################################

print('setup arguments from parser')

args = parser.parse_args()

#Define global variables.
lya = utils.lya_rest

if not args.pixels:
    args.pixels = list(range(12*args.nside**2))
N_pixels = len(args.pixels)

if args.save_suffix is not None:
    suffix = '_{}'.format(args.save_suffix)
else:
    suffix = ''

# TODO: print to confirm the arguments. e.g. "DLAs will be added"

if np.log2(args.nside)-int(np.log2(args.nside)) != 0:
    print('nside must be a power of 2!')
else:
    N_pix = 12*args.nside**2

#colours = ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'] * 2
colours = ['#F5793A','#A95AA1','#85C0F9','#0F2080'] * 2
fontsize = 16
figsize = (12, 6)
dpi = 80
#plt.rc('text', usetex=True)
plt.rc('font', size=fontsize)
plt.locator_params(axis='y', nbins=3)

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
    m = fits.open(args.base_dir+'/master.fits')
    if 'colorecell' in args.quantity:
        z = m['COSMO_COL'].data['Z']
        R = m['COSMO_COL'].data['R']
    else:
        z = m['COSMO_EXP'].data['Z']
        R = m['COSMO_EXP'].data['R']
    m.close()

    #Determine if we're looking at the Gaussian skewers.
    gaussian = ('colorecell' in args.quantity)

    dr_hMpc = (R[-1] - R[0])/(R.shape[0] - 1)

    #Function to get deltas and ivar from each pixel.
    def get_pixel_P1D(pixel,file_type='image'):
        if file_type == 'image':
            dirname = utils.get_dir_name(args.base_dir,pixel)
            filename = utils.get_file_name(dirname,'picca-'+args.quantity,args.nside,pixel,compressed=args.compressed_input)
            h = fits.open(filename)
            delta_rows = h[0].data.T
            ivar_rows = h[1].data.T
            z = 10**(h[2].data)/lya - 1
            h.close()
        elif file_type == 'delta':
            filename = args.base_dir + '/delta-{}.fits.gz'.format(pixel)
            h = fits.open(filename)
            for hdu in h[1:]:
                delta_rows = None
        Pk1D_results = {}
        for z_value in args.z_values:
            k, Pk, var = Pk1D.get_Pk1D(delta_rows,ivar_rows,dr_hMpc,z,z_value=z_value
                                    ,z_width=args.z_width,units=args.units,R1=args.smoothing_radius
                                    ,gaussian=gaussian)
            Pk1D_results[z_value] = {'k':k, 'Pk':Pk, 'var':var}
        return (pixel,Pk1D_results)

    tasks = [(pixel,args.file_type) for pixel in args.pixels]

    #Run the multiprocessing pool
    if __name__ == '__main__':
        pool = Pool(processes = args.nproc)
        results = []
        start_time = time.time()

        for task in tasks:
            pool.apply_async(get_pixel_P1D,task,callback=log_result,error_callback=log_error)

        pool.close()
        pool.join()

    #Combine the results from each pixel.
    for z_value in args.z_values:
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
        if z_value in args.z_values:
            Pk1D_results[float(hdu.name)] = {'k':hdu.data['k'], 'Pk':hdu.data['Pk'], 'var':hdu.data['var']}

################################################################################
"""
Save the data.
"""

#Save the data.
def save_P1D_values(Pk1D_results):

    header = fits.Header()
    header['units'] = args.units
    header['N_side'] = args.nside
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
    filename = 'Pk1D_data_{}_{}{}.fits'.format(args.quantity,N_pixels,suffix)
    hdulist.writeto(filename,overwrite=args.overwrite)
    hdulist.close

    return

if args.save_data:
    save_P1D_values(Pk1D_results)

################################################################################
"""
Make a plot.
"""

#Plot the power spectra
def plot_P1D_values(Pk1D_results,show_plot=True):
    fig = plt.figure(figsize=figsize, dpi= dpi, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')

    #Extract the z values and sort them.
    z_values = np.array([key for key in Pk1D_results.keys()])
    z_values.sort()

    #Set up to record the max/min values of P so we can set axis scales well.
    plot_min = 10**10.
    plot_max = 0.

    #For each z value, plot Pk1D and a model line.
    data_plots = {}
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

        #If desired plot points with error bars.
        if args.add_error_bars:

            #Choose the points to plot according to input number
            if args.N_k_values is not None:
                #Extract the k values we want to plot.
                N_rel = np.sum(k_rel)
                N_low = np.sum(k<args.k_min_plot)
                ind_orig = np.arange(N_low,N_low+N_rel,1)
                spacing = N_rel // args.N_k_values
                ind = ind_orig[(ind_orig % spacing) == 0]
                d = ax.errorbar(k[ind],to_plot_data[ind],yerr=to_plot_data_err[ind],label='z={}'.format(z_value),c=colour,fmt='o',zorder=3)
            else:
                d = ax.errorbar(k,to_plot_data,yerr=to_plot_data_err,label='z={}'.format(z_value),c=colour,fmt='o',zorder=3)

        #Otherwise, plot a line with error areas if required.
        elif args.add_error_areas:
            d = ax.plot(k,to_plot_data,label='z={}'.format(z_value),c=colour)
            ax.fill_between(k,to_plot_data+to_plot_data_err,to_plot_data-to_plot_data_err,color=colour,alpha=0.2,zorder=3)

        #Otherwise, just plot a line.
        else:
            d = ax.plot(k,to_plot_data,label='z={}'.format(z_value),c=colour,zorder=3)

        #Add the plotted data to the dictionary.
        data_plots[z_value] = d[0]

        #If we are plotting flux, add BOSS DR9 fitting function results.
        if 'flux' in args.quantity:
            model_Pk_kms = tuning.P1D_z_kms_PD2013(z_value,k)
            to_plot_model = model_Pk_kms * (k ** args.k_plot_power)
            plt.plot(k,to_plot_model,c=colour,linestyle='--',zorder=2)
            if args.add_model_interval_areas is not None:
                to_plot_lower = to_plot_model * (1 - args.add_model_interval_areas)
                to_plot_upper = to_plot_model * (1 + args.add_model_interval_areas)
                ax.fill_between(k,to_plot_upper,to_plot_lower,color=colour,alpha=0.2,zorder=1)
            elif args.add_model_interval_lines is not None:
                to_plot_lower = to_plot_model * (1 - args.add_model_interval_lines)
                to_plot_upper = to_plot_model * (1 + args.add_model_interval_lines)
                ax.plot(k,to_plot_upper,color=colour,linestyle=':',zorder=1)
                ax.plot(k,to_plot_lower,color=colour,linestyle=':',zorder=1)

        #Determine the minimum and maximum values plotted.
        #If we're plotting flux, use the model values for this.
        #Otherwise use the data values.
        if 'flux' in args.quantity:
            plot_min = np.minimum(plot_min,np.min(to_plot_model[k_rel]))
            plot_max = np.maximum(plot_max,np.max(to_plot_model[k_rel]))
        else:
            plot_min = np.minimum(plot_min,np.min(to_plot_data[k_rel]))
            plot_max = np.maximum(plot_max,np.max(to_plot_data[k_rel]))

    #Set the axis spacing correctly.
    space = 0.3
    ylim_lower = plot_min * (1 - space)
    ylim_upper = plot_max * (1 + space)
    ax.set_xlim(args.k_min_plot,args.k_max_plot)
    ax.set_ylim(ylim_lower,ylim_upper)

    #Add blank elements to describe data plot types in legend.
    descriptor_plots = {}
    if args.add_error_bars:
        d = ax.errorbar([0], [0], yerr=[0], fmt='o', color='gray', label=r'LyaColoRe')
    elif args.add_error_areas:
        d = ax.fill_between([0], [0], [0], color='gray', alpha=0.2, label=r'LyaColoRe')
    descriptor_plots['data'] = d

    #Add blank elements to describe model/fit/intervals plot types in legend.
    if 'flux' in args.quantity:
        err, = ax.plot([0], [0], color='gray', linestyle='--', label=r'BOSS DR9')
        if args.add_model_interval_lines is not None:
            err_i, = ax.plot([0], [0], color='gray', linestyle=':', label=r'DR9 $\pm$ 10%')
        elif args.add_model_interval_areas is not None:
            err_i = ax.fill_between([0], [0], [0], color='gray', alpha=0.2, label=r'DR9 $\pm$ 10%')
        descriptor_plots['error'] = err
        if args.add_model_interval_lines is not None or args.add_model_interval_areas is not None:
            descriptor_plots['error_interval'] = err_i

    #Add vertical line, and shade if desired.
    if args.add_vline is not None:
        ax.axvline(x=args.add_vline,linestyle='-.',color='darkgrey',zorder=4)
        if args.fill_right_vline:
            extend = 0.1
            k_shade = np.array([args.add_vline,args.k_max_plot*(1+extend)])
            y_shade_upper = np.array([ylim_upper*(1+extend)]*2)
            y_shade_lower = np.array([ylim_lower*(1-extend)]*2)
            #Line to fade the area with a translucent white overlay.
            #ax.fill_between(k_shade,y_shade_upper,y_shade_lower,linewidth=0.0,edgecolor=None,facecolor='white',alpha=0.5,zorder=4)
            #Line to fade the area with a translucent white overlay.
            #ax.fill_between(k_shade,y_shade_upper,y_shade_lower,linewidth=0.0,edgecolor=None,facecolor='grey',alpha=0.5,zorder=-1)
            #Line to hash the area with grey lines.
            ax.fill_between(k_shade,y_shade_upper,y_shade_lower,linewidth=0.0,edgecolor='darkgrey',facecolor='none',hatch='/',zorder=4)

    #Add a legend, ensuring everything is in the right order.
    legend_items = []
    legend_labels = []
    for z_value in np.sort(z_values)[::-1]:
        d = data_plots[z_value]
        legend_items += [d]
        legend_labels += [d.get_label()]
    for key in descriptor_plots.keys():
        legend_items += [descriptor_plots[key]]
        if type(descriptor_plots[key]) == list:
            legend_labels += [descriptor_plots[key][0].get_label()]
        else:
            legend_labels += [descriptor_plots[key].get_label()]
    ax.legend(legend_items,legend_labels,ncol=2,loc=0,facecolor='w')

    #Set yticks.
    if args.ysubs is not None:
        locator = matplotlib.ticker.LogLocator(subs=args.ysubs)
        ax.yaxis.set_major_locator(locator)
        thresh = np.ceil(np.log10(ylim_upper-ylim_lower))+1
        formatter = matplotlib.ticker.LogFormatterSciNotation(labelOnlyBase=False,minor_thresholds=(thresh,thresh))
        ax.yaxis.set_major_formatter(formatter)

    #Add appropriate labels.
    if args.k_plot_power == 1 :
        ylabel = r'$k\ P_{1\mathrm{D}}(k)$'
    elif args.k_plot_power == 2:
        ylabel = r'$k^{:d}\ P_{{1\mathrm{{D}}}}(k)\ [\mathrm{{s}}\ \mathrm{{km}}^{{-1}}]$'.format(int(args.k_plot_power))
    elif args.k_plot_power > 2 :
        ylabel = r'$k^{:d}\ P_{{1\mathrm{{D}}}}(k)\ [(\mathrm{{s}}\ \mathrm{{km}}^{{-1}})^{{{:d}}}]$'.format(int(args.k_plot_power),int(args.k_plot_power)-1,int(args.k_plot_power)-1)
    else:
        ylabel = r'$P_{1\mathrm{D}}\ [\mathrm{km\ s}^{{-1}}]$'
    ax.set_ylabel(ylabel,fontsize=fontsize)
    if args.units == 'km/s':
        ax.set_xlabel(r'$k\ [\mathrm{s\ km}^{-1}]$',fontsize=fontsize)
    elif args.units == 'Mpc/h':
        ax.set_xlabel(r'$k\ [h\ \mathrm{Mpc}^{-1}]$',fontsize=fontsize)
    filename = 'Pk1D_{}_{}.pdf'.format(args.quantity,N_pixels)
    plt.savefig('Pk1D_k{}_vs_k_{}{}.pdf'.format(int(args.k_plot_power),args.quantity,suffix))
    if show_plot:
        plt.show()
    return

plot_P1D_values(Pk1D_results,show_plot=args.show_plot)

################################################################################
"""
Celebrate!
"""
