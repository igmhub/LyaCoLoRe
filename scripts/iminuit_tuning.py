#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
import subprocess
from multiprocessing import Pool
import multiprocessing
import time
import glob
from iminuit import Minuit

from pyacolore import convert, Pk1D, utils, independent, tuning, simulation_data

lya = 1215.67

N_processes = 4
lambda_min = 3550.0
min_cat_z = 1.8
IVAR_cutoff = 1150.0

#Get the starting values of alpha, beta and sigma_G from file
#Decide which z values we are going to tune
z_values = [2.0]
z_width = 0.2

cell_size = 0.25 #Mpc/h

max_k = 0.005 #skm-1

#Open up the Gaussian colore files
base_file_location = '/Users/jfarr/Projects/test_data/test/'
#base_file_location = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/process_output_G_hZsmooth_4096_32_sr2.0_bm1_biasG18_picos_nside16_RSD'
N_side = 16

new_file_structure = '{}/{}/'               #pixel number//100, pixel number
new_filename_structure = '{}-{}-{}.fits'    #file type, nside, pixel number

input_format = 'gaussian_colore'

#get pixels from those directories created by make_master.py
dirs = glob.glob(base_file_location+new_file_structure.format('*','*'))
pixels = []
for dir in dirs[:2]:
    pixels += [int(dir[len(dir)-dir[-2::-1].find('/')-1:-1])]

# TODO: want to move this to tuning.py eventually
def measure_pixel_segment(pixel,z_value,alpha,beta,sigma_G_required,n,k1,A0):
    #print('start',pixel,z_value,alpha,beta,sigma_G_required,n,k1)

    location = base_file_location + '/' + new_file_structure.format(pixel//100,pixel)

    # TODO: this is a problem
    measured_SIGMA_G = 1.17

    #We work from the gaussian colore files made in 'pixelise gaussian skewers'.
    gaussian_filename = new_filename_structure.format('gaussian-colore',N_side,pixel)

    #Make a pixel object from it.
    data = simulation_data.SimulationData.get_gaussian_skewers_object(location+gaussian_filename,None,input_format,SIGMA_G=measured_SIGMA_G,IVAR_cutoff=IVAR_cutoff)

    #Determine the sigma_G to add
    extra_sigma_G = np.sqrt(sigma_G_required**2 - measured_SIGMA_G**2)

    #trim skewers
    data.trim_skewers(lambda_min,min_cat_z,extra_cells=1)

    #We extend the ranges of lambda a little to make sure RSDs are all accounted for.
    extra = 0.1
    lambda_min_val = lya*(1 + z_value - z_width*(1+extra)/2)
    lambda_max_val = lya*(1 + z_value + z_width*(1+extra)/2)
    data.trim_skewers(lambda_min_val,min_cat_z,lambda_max=lambda_max_val,whole_lambda_range=True)

    #add small scale fluctuations
    seed = int(str(N_side) + str(pixel))
    generator = np.random.RandomState(seed)
    data.add_small_scale_gaussian_fluctuations(cell_size,data.Z,np.ones(data.Z.shape[0])*extra_sigma_G,generator,amplitude=1.0,white_noise=False,lambda_min=0.0,IVAR_cutoff=lya,n=n,k1=k1,A0=A0) #n=0.7, k1=0.001 default

    #Convert to flux
    data.compute_physical_skewers()
    data.compute_tau_skewers(data.lya_absorber,alpha=np.ones(data.Z.shape[0])*alpha,beta=beta)
    data.add_RSDs(data.lya_absorber,np.ones(data.Z.shape[0])*alpha,beta,thermal=False)

    #Trim the skewers again to get rid of the additional cells
    lambda_min_val = lya*(1 + z_value - z_width/2)
    lambda_max_val = lya*(1 + z_value + z_width/2)
    data.trim_skewers(lambda_min_val,min_cat_z,lambda_max=lambda_max_val,whole_lambda_range=True)

    #Convert back to small cells for cf measurement
    #need a new function to merge cells back together

    ID = 0
    measurement = tuning.measurement(ID,z_value,z_width,data.N_qso,n,k1,alpha,beta,sigma_G_required,pixels=[pixel])

    measurement.add_mean_F_measurement(data)
    measurement.add_Pk1D_measurement(data)
    measurement.add_sigma_F_measurement(data)
    #print(measurement.sigma_F)

    measurement.add_mean_F_chi2(eps=0.1)
    measurement.add_Pk1D_chi2(max_k=max_k)
    measurement.add_sigma_F_chi2(eps=0.1)
    measurement.add_total_chi2()

    #print(tuning.get_flux_stats(sigma_G_required,alpha,beta,np.interp(z_value,data.Z,data.D)))

    return measurement

def f(alpha,beta,sigma_G,n,k1,A0):
    #beta=1.65
    ################################################################################

    """
    Define the multiprocessing tracking functions
    """

    #Define a progress-tracking function.
    def log_result(retval):

        results.append(retval)
        N_complete = len(results)
        N_tasks = len(tasks)

        #general.progress_bar(N_complete,N_tasks,start_time)
        return retval

    #Define an error-tracking function.
    def log_error(retval):
        print('Error:',retval)

    ################################################################################


    print('looking at params: a={:2.4f}, b={:2.4f}, sG={:2.4f}, n={:2.4f}, k1={:2.6f}, A0={:2.4f}'.format(alpha,beta,sigma_G,n,k1,A0))

    tasks = [(pixel,z_value,alpha,beta,sigma_G,n,k1,A0) for pixel in pixels for z_value in z_values]

    #Run the multiprocessing pool
    if __name__ == '__main__':
        pool = Pool(processes = N_processes)
        start_time = time.time()
        results = []

        for task in tasks:
            pool.apply_async(measure_pixel_segment,task,callback=log_result,error_callback=log_error)
            #results += [result]

        pool.close()
        pool.join()
    """

    results = []
    for task in tasks:
        results += [measure_pixel_segment(*task)]
    """

    measurement_set = tuning.measurement_set(measurements=results)
    combined_pixels_set = measurement_set.combine_pixels()

    Pk_kms_chi2 = 0.
    mean_F_chi2 = 0.
    overall_chi2 = 0.

    sigma_F_chi2 = 0.

    D = 0.3154712096325505 #at z=3
    predicted_flux_stats = tuning.get_flux_stats(sigma_G,alpha,beta,D)

    for m in combined_pixels_set.measurements:
        m.add_mean_F_chi2(eps=0.05)
        m.add_Pk1D_chi2(max_k=max_k)
        m.add_sigma_F_chi2(eps=0.05)
        m.add_total_chi2()
        Pk_kms_chi2 += m.Pk_kms_chi2
        mean_F_chi2 += m.mean_F_chi2
        overall_chi2 += m.total_chi2

        #print(' ')
        #print('mean F')
        #print('measured: {:2.4f}, model: {:2.4f}, predict: {:2.4f}'.format(m.mean_F,tuning.get_mean_F_model(m.z_value),predicted_flux_stats[0]))
        #print(' ')
        #print('sigma F')
        #print('measured: {:2.4f}, model: {:2.4f}, predict: {:2.4f}'.format(m.sigma_F,tuning.get_sigma_dF_P1D(m.z_value),predicted_flux_stats[1]))
        #print(' ')

        sigma_F_chi2 += m.sigma_F_chi2

    #chi2 = mean_F_chi2 + sigma_F_chi2
    #chi2 = mean_F_chi2 + Pk_kms_chi2
    chi2 = mean_F_chi2

    """
    combined_z_pixels_set = measurement_set.combine_zs()

    if len(combined_z_pixels_set.measurements) == 1:
        m = combined_z_pixels_set.measurements[0]
    """

    print('chi2: Pk {:2.4f}, mean F {:2.4f}, overall {:2.4f}'.format(Pk_kms_chi2,mean_F_chi2,overall_chi2))
    #print('chi2: sF {:2.4f}, mean F {:2.4f}, overall {:2.4f}'.format(sigma_F_chi2,mean_F_chi2,overall_chi2))
    print(' ')

    return chi2

t_kwargs = {'alpha' : 0.6095,    'error_alpha' : 0.05,   'limit_alpha' : (0., 10.),  'fix_alpha' : False,
            'beta' : 1.65,      'error_beta' : 0.05,    'limit_beta' : (0., 10.),   'fix_beta' : True,
            'sigma_G' : 5.237,  'error_sigma_G' : 0.05, 'limit_sigma_G' : (0., 20.),'fix_sigma_G' : False,
            }

t20_kwargs = {'alpha20' : 0.57,    'error_alpha20' : 0.05,   'limit_alpha20' : (0., 10.),  'fix_alpha20' : False,
            'beta20' : 1.65,      'error_beta20' : 0.05,    'limit_beta20' : (0., 10.),   'fix_beta20' : True,
            'sigma_G20' : 5.1,  'error_sigma_G20' : 0.05, 'limit_sigma_G20' : (0., 20.),'fix_sigma_G20' : False,
            }

t25_kwargs = {'alpha25' : 0.82,    'error_alpha25' : 0.05,   'limit_alpha25' : (0., 10.),  'fix_alpha25' : False,
            'beta25' : 1.65,      'error_beta25' : 0.05,    'limit_beta25' : (0., 10.),   'fix_beta25' : True,
            'sigma_G25' : 4.85,  'error_sigma_G25' : 0.05, 'limit_sigma_G25' : (0., 20.),'fix_sigma_G25' : False,
            }

t30_kwargs = {'alpha30' : 1.12,    'error_alpha30' : 0.05,   'limit_alpha30' : (0., 10.),  'fix_alpha30' : False,
            'beta30' : 1.65,      'error_beta30' : 0.05,    'limit_beta30' : (0., 10.),   'fix_beta30' : True,
            'sigma_G30' : 4.44,  'error_sigma_G30' : 0.05, 'limit_sigma_G30' : (0., 20.),'fix_sigma_G30' : False,
            }

s_kwargs = {'n'  : 1.051,       'error_n' : 0.05,       'limit_n' : (0., 10.),      'fix_n' : True, #0.9157
            'k1' : 0.005497,    'error_k1' : 0.00005,   'limit_k1' : (0., 0.1),     'fix_k1' : True, #0.003464
            'A0' : 58.6,        'error_A0' : 0.1,       'limit_A0' : (0., 200.),    'fix_A0' : True,
            }

minuit = Minuit(f,**t_kwargs,**s_kwargs)

minuit.print_param()
minuit.migrad() #ncall=20
minuit.print_param()

alpha = minuit.values['alpha']
beta = minuit.values['beta']
sigma_G = minuit.values['sigma_G']
n = minuit.values['n']
k1 = minuit.values['k1']
A0 = minuit.values['A0']

#Want to do a final run here
final_measurements = []
for pixel in pixels:
    for z_value in z_values:
        final_measurements += [measure_pixel_segment(pixel,z_value,alpha,beta,sigma_G,n,k1,A0)]
final_measurements = tuning.measurement_set(measurements=final_measurements)
final_measurements = final_measurements.combine_pixels()

#Plot a graph of mean F with redshift
mean_F = []
for m in final_measurements.measurements:
    mean_F += [m.mean_F]
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.plot(z_values,mean_F,label='mean flux',marker='o')
model_mean_F = tuning.get_mean_F_model(np.array(z_values))
plt.plot(z_values,model_mean_F,label='model',marker='o')
plt.fill_between(z_values,model_mean_F*1.1,model_mean_F*0.9,color=[0.5,0.5,0.5],alpha=0.5)
plt.grid()
plt.legend()
plt.show()

#Plot the power spectra
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
for m in final_measurements.measurements:
    plt.plot(m.k_kms,m.Pk_kms,label='z = {}'.format(m.z_value))
    model_Pk_kms = tuning.P1D_z_kms_PD2013(m.z_value,m.k_kms)
    plt.plot(m.k_kms,model_Pk_kms,label='model')
    plt.fill_between(m.k_kms,model_Pk_kms*1.1,model_Pk_kms*0.9,color=[0.5,0.5,0.5],alpha=0.5,label='model +/- 10%')
plt.axvline(x=max_k,color='k')
plt.semilogy()
plt.semilogx()
plt.grid()
plt.legend()
plt.savefig('Pk1D_z2.0.pdf')
plt.show()
