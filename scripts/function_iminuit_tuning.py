#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
from multiprocessing import Pool
import multiprocessing
import time
import glob
from iminuit import Minuit
import argparse

from lyacolore import convert, Pk1D, utils, independent, tuning, simulation_data

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Data options.
parser.add_argument('--base-dir', type = str, default = None, required=True,
                    help = 'Base directory for the input data')

parser.add_argument('--tuning-file-out', type = str, default = None, required=True,
                    help = 'Out file for the tuning data')

parser.add_argument('--plot-dir-out', type = str, default = None, required=False,
                    help = 'Out directory for the plots')

# Computational options.
parser.add_argument('--nproc', type = int, default = 1, required=False,
                    help = 'number of processes to use')

# Options for making skewers.
parser.add_argument('--pixels', type = int, default = None, required=False,
                    help = 'Which pixel numbers to use for input files', nargs='*')

parser.add_argument('--N-pixels', type = int, default = None, required=False,
                    help = 'Number of files to use as input')

parser.add_argument('--nside', type = int, default = 16, required=False,
                    help = 'HEALPix nside for output files (must be 2^n)')

parser.add_argument('--min-cat-z', type = float, default = 1.8, required=False,
                    help = 'minimum z of objects in catalog')

parser.add_argument('--seed', type = int, default = 16, required=False,
                    help = 'Random seed to use for the generation of random extra power.')

parser.add_argument('--cell-size', type = float, default = 0.25, required=False,
                    help = 'size in Mpc/h of output cells')

parser.add_argument('--lambda-min', type = float, default = 3550., required=False,
                    help = 'Minimum observed wavelength to use in the tuning process')

parser.add_argument('--lambda-rest-max', type = float, default = 1200., required=False,
                    help = 'Maximum rest-frame wavelength to use in the tuning process')

# Tuning options.
parser.add_argument('--z-values', type = float, default = [2.0,2.4,2.8,3.2], required=False,
                    help = 'which z values to measure at', nargs='*')

parser.add_argument('--z-width', type = float, default = 0.1, required=False,
                    help = 'Width of z bins')

# Output options.
parser.add_argument('--k-plot-max', type = float, default = 0.02, required=False,
                    help = 'max value of z to plot')

parser.add_argument('--overwrite', action="store_true", default = False, required=False,
                    help = 'overwrite existing files')

parser.add_argument('--compressed-input', action="store_true", default = False, required=False,
                    help = 'input files in format .fits.gz')

parser.add_argument('--show-plots', action="store_true", default = False, required=False,
                    help = 'show plots')

"""
# TODO: Implement these.
parser.add_argument('--start-from-file', type = str, default = None, required=False,
                    help = 'Tuning data file to use as a start point')

parser.add_argument('--start-from-random', action="store_true", default = False, required=False,
                    help = 'Tuning data file to use as a start point')

parser.add_argument('--nskewers', type = int, default = None, required=False,
                    help = 'number of skewers to process')

parser.add_argument('--downsampling', type = float, default = 1.0, required=False,
                    help = 'fraction by which to subsample the CoLoRe output')
"""
################################################################################

args = parser.parse_args()

lya = utils.lya_rest

#Admin options
save_plots = (args.plot_dir_out != None)
save_tuning = (args.tuning_file_out != None)

#Choose parameter values.
max_k = 0.01 #skm-1
eps_Pk1D = 0.1
eps_mean_F = 0.025
eps_bias_delta = 0.025
eps_bias_eta = 10**6#0.025
d_delta = 10.**-3
d_eta = 10**-9

# TODO: Enable tuning to be started from an ini file (needs a parser)
# TODO: Enable tuning to be started from an existing tuning file (with large errors)
# Not yet implemented
"""
#Get the starting values of the parameters from a tuning file.
default_tuning_file = '../input_files/tuning_data_with_bias_vel1.3_b1.65_lr1200.fits'
default_init_file = '../input_files/tuning.ini'
if args.start_from_file!=None and args.start_from_random:
    print('WARNING: Start file specified, and requested to start from random.')
    print(' -> Starting from {}'.format(args.tuning_file_start))
elif args.start_from_file==None and not args.start_from_random:
    print('WARNING: No start file specified, and not requested to start from random.')
    print(' -> Starting from {}'.format(default_tuning_file))
    args.tuning_file_start = default_tuning_file

if args.start_from_file!=None:
    iminuit_initial = tuning.iminuit_input_from_tuning_file(args.start_from_file)
elif args.start_from_random:
    iminuit_initial =
"""

# Choose tuning parameter initial values.
initial_C0 = 1.482229863221668
initial_C1 = 4.5
initial_C2 = 0.0
initial_texp = 1.65
initial_D0 = 6.018308640829534
initial_D1 = 0.2756162052010332
initial_D2 = 0.0
initial_n = 0.7318824370864454
initial_k1 = 0.0341049400675243
initial_R = 25.0 #kms-1
initial_a_v = 1.3

# Choose parameters to fix.
fix_all = True
fix_C0 = False
fix_C1 = True
fix_C2 = True
fix_texp = True
fix_D0 = False
fix_D1 = False
fix_D2 = True
fix_n = False
fix_k1 = False
fix_R = True
fix_a_v = True

# Make dictionaries for input to iminuit.
tau0_kwargs = {'C0' : initial_C0,      'error_C0' : 1.0,   'fix_C0' : fix_all|fix_C0,     'limit_C0' : (0., 100.),
               'C1' : initial_C1,      'error_C1' : 1.0,   'fix_C1' : fix_all|fix_C1,     #'limit_C1' : (0., 20.),
               'C2' : initial_C2,      'error_C2' : 1.0,   'fix_C2' : fix_all|fix_C2,     #'limit_C2' : (0., 20.),
               }

texp_kwargs = {'texp' : initial_texp,  'error_texp' : 1.0, 'fix_texp' : fix_all|fix_texp, 'limit_texp' : (0.,5.)
               }

seps_kwargs = {'D0' : initial_D0,     'error_D0' : 1.0,   'fix_D0' : fix_all|fix_D0,     'limit_D0' : (0., 100.),
               'D1' : initial_D1,     'error_D1' : 0.2,   'fix_D1' : fix_all|fix_D1,     #'limit_D1' : (0., 20.),
               'D2' : initial_D2,     'error_D2' : 1.0,   'fix_D2' : fix_all|fix_D2,     #'limit_D2' : (0., 20.),
               }

s_kwargs = {'n'  : initial_n,       'error_n' : 1.0,    'fix_n' : fix_all|fix_n,       'limit_n' : (-2., 10.),
            'k1' : initial_k1,      'error_k1' : 0.001, 'fix_k1' : fix_all|fix_k1,     'limit_k1' : (0., 0.1),
            }

other_kwargs = {'R'  : initial_R,    'error_R' : 1.0,   'fix_R' : fix_all|fix_R,       'limit_R' : (0., 1000.),
                'a_v': initial_a_v,  'error_a_v' : 0.1, 'fix_a_v' : fix_all|fix_a_v,   'limit_a_v' : (0., 2.0),
                'return_measurements'  : False,    'fix_return_measurements' : True,
                'errordef'             : 1,
                }

#Colours for the plot
colours = ['C0','C1','C2','C3','C4','C5','C6']

#Function to produce a parameter as a function of z.
def get_parameter(z,A0,A1,A2,z0=3.0):
    x = (1 + z)/(1 + z0)
    alpha = np.exp(utils.quadratic_log(x,A0,A1,A2))
    return alpha

#Check that the output tuning file has the right extension.
if (args.tuning_file_out[-8:] != '.fits.gz') and (args.tuning_file_out[-5:] != '.fits'):
    raise NameError('Output tuning file is not .fits or .fits.gz')

#Get the location to save the plots if none is given.
if args.plot_dir_out is None:
    print('No plot directory specified: plots will not be saved.')

#Get the list of input pixels if one is not given.
if (args.pixels is None) and (args.N_pixels is None):
    print('ALERT: Neither input pixels nor number of input pixels specified: using pixels 0 to 63.')
    args.pixels = list(range(64))
elif (args.pixels is None):
    print('ALERT: Input pixels not specified: using pixels 0 to {}.'.format(args.N_pixels))
    args.pixels = list(range(args.N_pixels))

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
    return retval

#Define an error-tracking function.
def log_error(retval):
    print('Error:',retval)

################################################################################

# TODO: want to move this to tuning.py eventually
def measure_pixel_segment(pixel,C0,C1,C2,texp,D0,D1,D2,n,k1,R_kms,a_v,RSD_weights,prep=False):
    t = time.time()
    seed = int(pixel * 10**5 + args.seed)

    #print('start pixel {} at {}'.format(pixel,time.ctime()))

    #Get the filename of the gaussian skewer.
    location = utils.get_dir_name(args.base_dir,pixel)
    filename = utils.get_file_name(location,'gaussian-colore',args.nside,pixel,compressed=args.compressed_input)

    #Make a pixel object from it.
    data = simulation_data.SimulationData.get_gaussian_skewers_object(filename,None,'gaussian_colore',IVAR_cutoff=args.lambda_rest_max)
    #print('{:3.2f} checkpoint sim_dat'.format(time.time()-t))
    t = time.time()

    #Scale the RSD skewers.
    data.scale_velocities(a_v=a_v)

    #Get the transformation for the current set of input parameters.
    transformation = tuning.Transformation()
    def f_tau0_z(z):
        return get_parameter(z,C0,C1,C2)
    def f_texp_z(z):
        return get_parameter(z,texp,0.,0.)
    def f_seps_z(z):
        return get_parameter(z,D0,D1,D2)
    transformation.add_parameters_from_functions(f_tau0_z,f_texp_z,f_seps_z)
    data.add_transformation(transformation)

    #trim skewers to the minimal length
    lambda_buffer = 100. #Angstroms
    z_lower_cut = np.min(args.z_values) - args.z_width/2.
    z_upper_cut = np.max(args.z_values) + args.z_width/2.
    lambda_min_val = np.min([args.lambda_min,lya*(1 + z_lower_cut)]) - lambda_buffer
    lambda_max_val = lya*(1 + z_upper_cut) + lambda_buffer
    data.trim_skewers(lambda_min_val,args.min_cat_z,lambda_max=lambda_max_val,whole_lambda_range=False)

    #Add small scale fluctuations to the skewers.
    generator = np.random.RandomState(seed)
    data.add_small_scale_gaussian_fluctuations(args.cell_size,generator,white_noise=False,lambda_min=0.0,IVAR_cutoff=args.lambda_rest_max,n=n,k1=k1,R_kms=R_kms)

    #print('{:3.2f} checkpoint extra power'.format(time.time()-t))
    t = time.time()

    #Copmute the physical skewers
    data.compute_physical_skewers()

    #Compute the tau skewers and add RSDs
    data.compute_tau_skewers(data.lya_absorber)
    #print('{:3.2f} checkpoint tau'.format(time.time()-t))
    t = time.time()

    if prep:
        data.compute_RSD_weights(thermal=False)

        #print(pixel,'{:3.2f} checkpoint RSD weights measured'.format(time.time()-t))
        t = time.time()

        #b_eta_weights_dict = data.get_bias_eta_RSD_weights(args.z_values,d=d_eta,z_width=args.z_width,lambda_buffer=lambda_buffer)
        b_eta_weights_dict = None

        #print(pixel,'{:3.2f} checkpoint b_eta weights measured'.format(time.time()-t))
        t = time.time()

        return (pixel,data.RSD_weights,b_eta_weights_dict)
    else:
        RSD_weights = RSD_weights_dict[pixel]
        bias_eta_weights = bias_eta_weights_dict[pixel]
        data.add_all_RSDs(thermal=False,weights=RSD_weights)

        #print('{:3.2f} checkpoint RSDs'.format(time.time()-t))
        t = time.time()

        measurements = []
        times_m = np.zeros(6)
        for z_value in args.z_values:
            ID = n
            t_m = time.time()
            measurement = tuning.function_measurement(ID,z_value,args.z_width,data.N_qso,n,k1,C0,C1,C2,texp,D0,D1,D2,pixels=[pixel])
            times_m[0] += time.time() - t_m
            t_m = time.time()
            measurement.add_mean_F_measurement(data)
            times_m[1] += time.time() - t_m
            t_m = time.time()
            measurement.add_Pk1D_measurement(data)
            times_m[2] += time.time() - t_m
            t_m = time.time()
            #measurement.add_sigma_dF_measurement(data)
            times_m[3] += time.time() - t_m
            t_m = time.time()
            measurement.add_bias_delta_measurement(data,d=d_delta,weights=RSD_weights)
            times_m[4] += time.time() - t_m
            t_m = time.time()
            #measurement.add_bias_eta_measurement(data,d=d_eta,weights_dict=bias_eta_weights,lambda_buffer=lambda_buffer)
            times_m[5] += time.time() - t_m
            measurements += [measurement]

        #print('{:3.2f} checkpoint measurements'.format(time.time()-t))
        #print('--> measurement_times: {:3.2f}, {:3.2f}, {:3.2f}, {:3.2f}, {:3.2f}, {:3.2f}'.format(times_m[0],times_m[1],times_m[2],times_m[3],times_m[4],times_m[5]))

        return measurements

################################################################################
"""
Produce the RSD weights matrices
"""

print('Producing preparatory RSD maps:')
tasks = [(pixel,initial_C0,initial_C1,initial_C2,initial_texp,initial_D0,initial_D1,initial_D2,initial_n,initial_k1,initial_R,initial_a_v,None,True) for pixel in pixels]

if __name__ == '__main__':
    pool = Pool(processes = args.nproc)
    start_time = time.time()
    results = []

    for task in tasks:
        pool.apply_async(measure_pixel_segment,task,callback=log_result,error_callback=log_error)

    pool.close()
    pool.join()

RSD_weights_dict = {}
bias_eta_weights_dict = {}
for result in results:
    RSD_weights_dict[result[0]] = result[1]
    bias_eta_weights_dict[result[0]] = result[2]
print('done!\n')

################################################################################

#Define the function which we seek to minimise.
def f(C0,C1,C2,texp,D0,D1,D2,n,k1,R,a_v,return_measurements=False):

    #Re-define the multiprocessing tracking functions without the progress bar
    def log_result(retval):
        results.append(retval)
        return retval
    def log_error(retval):
        print('Error:',retval)

    print('starting at',time.ctime())
    print('looking at params: C=({:2.4f},{:2.4f},{:2.4f}), texp={:1.2f},  D=({:2.4f},{:2.4f},{:2.4f}), n={:2.4f}, k1={:2.6f}'.format(C0,C1,C2,texp,D0,D1,D2,n,k1))

    tasks = [(pixel,C0,C1,C2,texp,D0,D1,D2,n,k1,R,a_v,None) for pixel in pixels]

    #Run the multiprocessing pool
    if __name__ == '__main__':
        pool = Pool(processes = args.nproc)
        start_time = time.time()
        results = []

        for task in tasks:
            pool.apply_async(measure_pixel_segment,task,callback=log_result,error_callback=log_error)

        pool.close()
        pool.join()

    measurements = []
    for result in results:
        measurements += result
    measurement_set = tuning.measurement_set(measurements=measurements)
    combined_pixels_set = measurement_set.combine_pixels()

    Pk_kms_chi2 = 0.
    mean_F_chi2 = 0.
    bias_delta_chi2 = 0.
    bias_eta_chi2 = 0.

    for m in combined_pixels_set.measurements:
        m.add_mean_F_chi2(eps=eps_mean_F)
        m.add_Pk1D_chi2(max_k=max_k,denom="npower_cutoff",eps=eps_Pk1D)
        m.add_bias_delta_chi2(eps=eps_bias_delta)
        #m.add_bias_eta_chi2(eps=eps_bias_eta)

        Pk_kms_chi2 += m.Pk_kms_chi2
        mean_F_chi2 += m.mean_F_chi2
        bias_delta_chi2 += m.bias_delta_chi2
        #bias_eta_chi2 += m.bias_eta_chi2
        #print('z =',m.z_value,'number of k values =',m.k_kms.shape,'number of k values < max k =',np.sum(m.k_kms<max_k))

    chi2 = mean_F_chi2 + Pk_kms_chi2 + bias_delta_chi2 + bias_eta_chi2
    log_text = 'chi2: Pk {:2.4f}, mean F {:2.4f}, b_del {:2.4f}, b_eta {:2.4f}, overall {:2.4f}'.format(Pk_kms_chi2,mean_F_chi2,bias_delta_chi2,bias_eta_chi2,chi2)

    print(log_text)
    print(' ')

    with open("parameter_log.txt","a") as f:
        txt = str(time.ctime()+'\n')
        txt += 'C0:{}, C1:{}, C2:{}, D0:{}, D1:{}, D2:{}, n:{}, k1:{}, texp:{}\n'.format(C0,C1,C2,D0,D1,D2,n,k1,texp)
        txt += log_text + '\n\n'
        f.write(txt)
        f.close()
        best = chi2

    if return_measurements:
        return chi2, combined_pixels_set
    else:
        return chi2

minuit = Minuit(f,**tau0_kwargs,**texp_kwargs,**seps_kwargs,**s_kwargs,**other_kwargs)

if not fix_all:
    minuit.print_param()
    minuit.migrad() #ncall=20
    minuit.print_param()

C0 = minuit.values['C0']
C1 = minuit.values['C1']
C2 = minuit.values['C2']
texp = minuit.values['texp']
D0 = minuit.values['D0']
D1 = minuit.values['D1']
D2 = minuit.values['D2']
n = minuit.values['n']
k1 = minuit.values['k1']
R = minuit.values['R']
a_v = minuit.values['a_v']

print(minuit.values)

#Function to define the output tuning file.
def save_tuning_file(filename,overwrite=False):
    z = np.linspace(1.6,4.0,2401)
    alpha_arr = get_parameter(z,C0,C1,C2)
    texp_arr = np.ones(alpha_arr.shape)*texp
    sigma_G_arr = get_parameter(z,D0,D1,D2)

    header = fits.Header()
    header['C0'] = C0
    header['C1'] = C1
    header['C2'] = C2
    header['D0'] = D0
    header['D1'] = D1
    header['D2'] = D2
    header['n'] = n
    header['k1'] = k1
    header['R'] = R
    header['a_v'] = a_v

    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)

    dtype = [('z', 'f4'), ('alpha', 'f4'), ('texp', 'f4'), ('sigma_G', 'f4')]
    data = np.array(list(zip(z,alpha_arr,texp_arr,sigma_G_arr)),dtype=dtype)
    hdu_tuning = fits.BinTableHDU(data,header=header,name='TUNING')

    hdulist = fits.HDUList([prihdu,hdu_tuning])
    hdulist.writeto(filename,overwrite=overwrite)
    hdulist.close

    return

if save_tuning:
    save_tuning_file(args.tuning_file_out,overwrite=args.overwrite)

#Do a final run to get the measurements.
final_chi2,final_measurements = f(C0,C1,C2,texp,D0,D1,D2,n,k1,R,a_v,return_measurements=True)

#Plot a graph of mean_F, bias_delta or bias_eta with redshift
def plot_scalar_measurement_values(m_set,data_type='mean_F',show_plot=True,save_plot=False):
    data = []
    for m in m_set.measurements:
        if data_type == 'mean_F':
            data += [(m.z_value,m.mean_F)]
        elif data_type == 'bias_delta':
            data += [(m.z_value,m.bias_delta)]
        elif data_type == 'bias_eta':
            data += [(m.z_value,m.bias_eta)]
    dtype = [('z', 'd'), ('data', 'd')]
    data = np.array(data,dtype=dtype)
    data = np.sort(data,order=['z'])

    if data_type == 'mean_F':
        text = r'$\bar{F}$'
        model = tuning.get_mean_F_model(np.array(args.z_values))
    elif data_type == 'bias_delta':
        text = r'$b_\delta$'
        model = tuning.get_bias_delta_model(np.array(args.z_values))
    elif data_type == 'bias_eta':
        text = r'$b_\eta$'
        model = tuning.get_bias_eta_model(np.array(args.z_values))

    plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
    plt.plot(data['z'],data['data'],label=text+' mocks',marker='o')

    plt.plot(args.z_values,model,label=text+' model',marker='x',c=[0.5,0.5,0.5])
    plt.fill_between(args.z_values,model*1.1,model*0.9,color=[0.5,0.5,0.5],alpha=0.3)
    plt.grid()
    plt.legend()

    z = data['z']
    min_z = np.min(z) - (np.max(z) - np.min(z))*0.1
    max_z = np.max(z) + (np.max(z) - np.min(z))*0.1
    plt.xlim(min_z,max_z)
    plt.xlabel(r'$z$',fontsize=12)
    plt.ylabel(text,fontsize=12)
    if save_plot:
        plt.savefig(args.plot_dir_out + data_type + '.pdf')
    if show_plot:
        plt.show()
    return

#Plot the power spectra
def plot_P1D_values(m_set,show_plot=True,save_plot=False):
    plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
    max_power_plot = 0.
    min_power_plot = 10.**10
    for m in m_set.measurements:
        colour = colours[np.searchsorted(args.z_values,m.z_value)]
        plt.plot(m.k_kms,m.Pk_kms,label=r'$z={}$'.format(m.z_value),c=colour)
        model_Pk_kms = tuning.P1D_z_kms_PD2013(m.z_value,m.k_kms)
        plt.plot(m.k_kms,model_Pk_kms,c=colour,linestyle=':')#,label='z={} DR9'.format(m.z_value))
        m.add_Pk1D_chi2(max_k=max_k,denom="npower_cutoff")
        eps = m.Pk_kms_chi2_eps
        #plt.plot(m.k_kms,model_Pk_kms*0.9,color=[0.5,0.5,0.5],alpha=0.5)
        #plt.plot(m.k_kms,model_Pk_kms*1.1,color=[0.5,0.5,0.5],alpha=0.5)
        #plt.fill_between(m.k_kms,model_Pk_kms*1.1,model_Pk_kms*0.9,color=[0.5,0.5,0.5],alpha=0.3)#,label='DR9 +/- 10%')
        lower = np.maximum(np.ones_like(model_Pk_kms)*10**(-6),model_Pk_kms * (1. - eps))
        upper = model_Pk_kms * (1. + eps)
        #plt.plot(m.k_kms,upper,c='k',linestyle='dashed')
        #plt.plot(m.k_kms,lower,c='k',linestyle='dashed')
        plt.fill_between(m.k_kms,upper,lower,color=[0.8,0.8,0.8],alpha=0.3)
        k_valid = m.k_kms<args.k_plot_max
        max_power_plot = np.max((max_power_plot,np.max(model_Pk_kms[k_valid])))
        min_power_plot = np.min((min_power_plot,np.min(model_Pk_kms[k_valid])))
    #plt.title('C=({:2.4f},{:2.4f},{:2.4f}), D=({:2.4f},{:2.4f},{:2.4f}), n={:2.4f}, k1={:2.6f}'.format(m.C0,m.C1,m.C2,m.D0,m.D1,m.D2,m.n,m.k1),fontsize=12)
    #plt.axvline(x=max_k,color='k')
    plt.semilogy()
    plt.semilogx()
    ylim_lower = min_power_plot * 0.8
    ylim_upper = max_power_plot * 1.2
    plt.ylim(ylim_lower,ylim_upper)
    plt.xlim(xmax=args.k_plot_max)
    plt.grid()
    plt.legend(fontsize=12)
    plt.ylabel(r'$P_{1D}$',fontsize=12)
    plt.xlabel(r'$k\ /\ (kms^{-1})^{-1}$',fontsize=12)
    if save_plot:
        plt.savefig(args.plot_dir_out + '/Pk1D.pdf')
    if show_plot:
        plt.show()
    return

#Plot parameter values
def plot_parameter_values(minuit_results,z_min=1.8,z_max=4.0,z_size=0.01,show_plot=True,save_plot=False):
    C0 = minuit.values['C0']
    C1 = minuit.values['C1']
    C2 = minuit.values['C2']
    D0 = minuit.values['D0']
    D1 = minuit.values['D1']
    D2 = minuit.values['D2']
    n = minuit.values['n']
    k1 = minuit.values['k1']

    z = np.linspace(z_min,z_max,(z_max-z_min)/z_size+1)
    alpha = get_parameter(z,C0,C1,C2)
    sigma_G = get_parameter(z,D0,D1,D2)

    """
    #HACK
    manual_cells = z<2.5
    transfer_cell = np.searchsorted(z,2.5)
    grad = (alpha[transfer_cell] - alpha[transfer_cell + 1])/(z[transfer_cell + 1] - z[transfer_cell])
    manual_alphas = alpha[transfer_cell] + grad*(z[transfer_cell] - z)
    alpha = alpha - manual_cells * (alpha - manual_alphas)
    """

    plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
    plt.plot(z,alpha,label=r'$\alpha$')
    plt.plot(z,sigma_G,label=r'$\sigma_\epsilon$')
    plt.title('C=({:2.4f},{:2.4f},{:2.4f}), D=({:2.4f},{:2.4f},{:2.4f}), n={:2.4f}, k1={:2.6f}'.format(C0,C1,C2,D0,D1,D2,n,k1),fontsize=12)
    plt.grid()
    plt.legend(fontsize=12)
    plt.ylabel('Parameters',fontsize=12)
    plt.xlabel(r'$z$',fontsize=12)
    min_z = np.min(z) - (np.max(z) - np.min(z))*0.1
    max_z = np.max(z) + (np.max(z) - np.min(z))*0.1
    plt.xlim(min_z,max_z)
    if save_plot:
        plt.savefig(args.plot_dir_out + '/parameters.pdf')
    if show_plot:
        plt.show()
    return

#Make the relevant plots.
plot_scalar_measurement_values(final_measurements,data_type='mean_F',show_plot=args.show_plots,save_plot=save_plots)
plot_scalar_measurement_values(final_measurements,data_type='bias_delta',show_plot=args.show_plots,save_plot=save_plots)
plot_scalar_measurement_values(final_measurements,data_type='bias_eta',show_plot=args.show_plots,save_plot=save_plots)
plot_P1D_values(final_measurements,show_plot=args.show_plots,save_plot=save_plots)
plot_parameter_values(minuit,z_min=min(args.z_values)-args.z_width/2.,z_max=max(args.z_values)+args.z_width/2.,show_plot=args.show_plots,save_plot=save_plots)
