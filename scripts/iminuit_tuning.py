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

from pyacolore import convert, Pk1D, utils, independent, tuning, simulation_data

lya = 1215.67

N_processes = 32
lambda_min = 3550.0
min_cat_z = 1.8
IVAR_cutoff = 1150.0

#Get the starting values of alpha, beta and sigma_G from file
#Decide which z values we are going to tune
z_value = 3.0
z_width = 0.2

cell_size = 0.25 #Mpc/h

max_k = 0.005 #skm-1

#Open up the Gaussian colore files
base_file_location = '/Users/jfarr/Projects/test_data/test/'
base_file_location = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/process_output_G_hZsmooth_4096_32_sr2.0_bm1_biasG18_picos_newNz_mpz0_seed1003_123_nside16/'
N_side = 16

new_file_structure = '{}/{}/'               #pixel number//100, pixel number
new_filename_structure = '{}-{}-{}.fits'    #file type, nside, pixel number

input_format = 'gaussian_colore'

#get pixels from those directories created by make_master.py
dirs = glob.glob(base_file_location+new_file_structure.format('*','*'))
pixels = []
for dir in dirs[:32]:
    ending = dir[len(dir)-dir[-2::-1].find('/')-1:-1]
    if ending != 'logs':
        pixels += [int(ending)]

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

    #print('before SSP')
    mean_G = np.average(data.GAUSSIAN_DELTA_rows)
    #print('gaussian field mean is',mean_G)
    mean_G2 = np.average(data.GAUSSIAN_DELTA_rows**2)
    sigma_G_simple = np.sqrt(mean_G2 - mean_G**2)
    #print('gaussian field simple sigma is',sigma_G_simple)
    #print('gaussian field std is',np.std(data.GAUSSIAN_DELTA_rows))
    dkms_dhMpc = utils.get_dkms_dhMpc(z_value)
    dv_kms = dkms_dhMpc*(data.R[-1] - data.R[0])/(data.N_cells - 1)
    k = np.fft.rfftfreq(data.N_cells)*2*np.pi/dv_kms
    pk_rows = np.fft.rfft(data.GAUSSIAN_DELTA_rows,axis=1) / np.sqrt(data.N_cells/dv_kms)
    pk_rows = np.abs(pk_rows)**2
    pk_measured = np.average(pk_rows,axis=0)
    #print('gaussian field Pk sigma is',np.sqrt(np.trapz(pk_measured,k)/(np.pi)))
    #print(' ')

    #add small scale fluctuations
    seed = int(str(N_side) + str(pixel))
    generator = np.random.RandomState(seed)
    data.add_small_scale_gaussian_fluctuations(cell_size,data.Z,np.ones(data.Z.shape[0])*extra_sigma_G,generator,amplitude=1.0,white_noise=False,lambda_min=0.0,IVAR_cutoff=lya,n=n,k1=k1,A0=A0) #n=0.7, k1=0.001 default

    #print('after SSP')
    mean_G = np.average(data.GAUSSIAN_DELTA_rows)
    #print('gaussian field mean is',mean_G)
    mean_G2 = np.average(data.GAUSSIAN_DELTA_rows**2)
    sigma_G_simple = np.sqrt(mean_G2 - mean_G**2)
    #print('gaussian field simple sigma is',sigma_G_simple)
    #print(' ')

    #Convert to flux
    data.compute_physical_skewers()

    #print('physical')
    mean_D = np.average(data.DENSITY_DELTA_rows)
    #print('density delta field mean is',mean_D)
    mean_D2 = np.average(data.GAUSSIAN_DELTA_rows**2)
    sigma_D_simple = np.sqrt(mean_D2 - mean_D**2)
    #print('density delta field simple sigma is',sigma_D_simple)
    int_lim = sigma_G_simple*10.
    delta_G_integral = np.linspace(-int_lim,int_lim,10**4)
    delta_G_integral = np.reshape(delta_G_integral,(1,delta_G_integral.shape[0]))
    prob_delta_G = (1/((np.sqrt(2*np.pi))*sigma_G_simple))*np.exp(-(delta_G_integral**2)/(2*(sigma_G_simple**2)))
    D = np.interp(z_value,data.Z,data.D) #* np.ones_like(delta_G_integral)
    density_delta_integral = convert.gaussian_to_lognormal_delta(delta_G_integral,sigma_G_simple,D)
    mean_D = np.trapz(prob_delta_G * density_delta_integral,delta_G_integral)[0]
    sigma_D = np.sqrt(np.trapz(prob_delta_G * (density_delta_integral)**2,delta_G_integral)[0])
    #print('predicted density delta mean is',mean_D)
    #print('predicted density delta sigma is',sigma_D)
    #print(' ')

    data.compute_tau_skewers(data.lya_absorber,alpha=np.ones(data.Z.shape[0])*alpha,beta=beta)
    data.add_RSDs(data.lya_absorber,np.ones(data.Z.shape[0])*alpha,beta,thermal=False)

    #Trim the skewers again to get rid of the additional cells
    lambda_min_val = lya*(1 + z_value - z_width/2)
    lambda_max_val = lya*(1 + z_value + z_width/2)
    data.trim_skewers(lambda_min_val,min_cat_z,lambda_max=lambda_max_val,whole_lambda_range=True)

    #Convert back to small cells for cf measurement
    #need a new function to merge cells back together

    ID = n
    measurement = tuning.measurement(ID,z_value,z_width,data.N_qso,n,k1,alpha,beta,sigma_G_required,pixels=[pixel])

    measurement.add_mean_F_measurement(data)
    measurement.add_Pk1D_measurement(data)
    measurement.add_sigma_F_measurement(data)
    #print(measurement.sigma_F)

    measurement.add_mean_F_chi2(eps=0.05)
    measurement.add_Pk1D_chi2(max_k=max_k)
    measurement.add_sigma_F_chi2(eps=0.1)
    measurement.add_total_chi2()

    #print(tuning.get_flux_stats(sigma_G_required,alpha,beta,np.interp(z_value,data.Z,data.D)))


    return measurement

def f(alpha,beta,sigma_G,n,k1,A0,return_measurements=False):

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

    tasks = [(pixel,z_value,alpha,beta,sigma_G,n,k1,A0) for pixel in pixels]

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

    measurement_set = tuning.measurement_set(measurements=results)
    combined_pixels_set = measurement_set.combine_pixels()

    Pk_kms_chi2 = 0.
    mean_F_chi2 = 0.
    overall_chi2 = 0.

    sigma_F_chi2 = 0.

    D = 0.4172239560422373 #at z=2.0
    #D = 0.35947529665680117 #at z=2.5
    #D = 0.3154712096325505 #at z=3

    predicted_flux_stats = tuning.get_flux_stats(sigma_G,alpha,beta,D)

    for m in combined_pixels_set.measurements:
        m.add_mean_F_chi2(eps=0.05)
        m.add_Pk1D_chi2(max_k=max_k)#,denom="npower")
        m.add_sigma_F_chi2(eps=0.05)
        m.add_total_chi2()
        Pk_kms_chi2 += m.Pk_kms_chi2
        mean_F_chi2 += m.mean_F_chi2
        overall_chi2 += m.total_chi2

        print(' ')
        print('mean F')
        print('measured: {:2.4f}, model: {:2.4f}, predict: {:2.4f}'.format(m.mean_F,tuning.get_mean_F_model(m.z_value),predicted_flux_stats[0]))
        print(' ')
        print('sigma F')
        print('measured: {:2.4f}, model: {:2.4f}, predict: {:2.4f}'.format(m.sigma_F,tuning.get_sigma_dF_P1D(m.z_value),predicted_flux_stats[1]))
        print(' ')

        sigma_F_chi2 += m.sigma_F_chi2

    #chi2 = mean_F_chi2 + sigma_F_chi2
    chi2 = mean_F_chi2 + Pk_kms_chi2
    #chi2 = mean_F_chi2

    """
    combined_z_pixels_set = measurement_set.combine_zs()

    if len(combined_z_pixels_set.measurements) == 1:
        m = combined_z_pixels_set.measurements[0]
    """

    print('chi2: Pk {:2.4f}, mean F {:2.4f}, overall {:2.4f}'.format(Pk_kms_chi2,mean_F_chi2,overall_chi2))
    #print('chi2: sF {:2.4f}, mean F {:2.4f}, overall {:2.4f}'.format(sigma_F_chi2,mean_F_chi2,overall_chi2))
    print(' ')

    if return_measurements:
        return chi2, combined_pixels_set
    else:
        return chi2

t_kwargs = {'alpha' : 2.30,    'error_alpha' : 0.05,   'limit_alpha' : (0., 30.),  'fix_alpha' : True,
            'beta' : 1.65,      'error_beta' : 0.05,    'limit_beta' : (0., 20.),   'fix_beta' : True,
            'sigma_G' : 6.86,  'error_sigma_G' : 0.05, 'limit_sigma_G' : (0., 20.),'fix_sigma_G' : True,
            }

s_kwargs = {'n'  : 1.52,       'error_n' : 0.05,       'limit_n' : (0., 10.),      'fix_n' : True,
            'k1' : 0.0166,    'error_k1' : 0.0005,   'limit_k1' : (0., 0.1),     'fix_k1' : True,
            'A0' : 58.6,        'error_A0' : 0.1,       'limit_A0' : (0., 200.),    'fix_A0' : True,
            }

other_kwargs = {'return_measurements'  : False,    'fix_return_measurements' : True,
            }

minuit = Minuit(f,**t_kwargs,**s_kwargs,**other_kwargs)

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
#final_measurements = []
final_chi2,final_measurements = f(alpha=alpha,beta=beta,sigma_G=sigma_G,n=n,k1=k1,A0=A0,return_measurements=True)
#for pixel in pixels:
#    final_measurements += [measure_pixel_segment(pixel,z_value,alpha,beta,sigma_G,n,k1,A0)]
#final_measurements = tuning.measurement_set(measurements=final_measurements)
#final_measurements = final_measurements.combine_pixels()


"""
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
"""

#Plot the power spectra
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
for m in final_measurements.measurements:
    plt.plot(m.k_kms,m.Pk_kms,label='z = {}'.format(m.z_value))
    model_Pk_kms = tuning.P1D_z_kms_PD2013(m.z_value,m.k_kms)
    plt.plot(m.k_kms,model_Pk_kms,label='DR9 fitting function')
    m.add_Pk1D_chi2(max_k=max_k)#,denom="npower")
    eps = m.Pk_kms_chi2_eps
    plt.fill_between(m.k_kms,model_Pk_kms*1.1,model_Pk_kms*0.9,color=[0.5,0.5,0.5],alpha=0.5,label='DR9 +/- 10%')
    lower = np.maximum(np.ones_like(model_Pk_kms)*10**(-6),model_Pk_kms * (1. - eps))
    upper = model_Pk_kms * (1. + eps)
    plt.plot(m.k_kms,upper,c='k',linestyle='dashed')
    plt.plot(m.k_kms,lower,c='k',linestyle='dashed')
    plt.fill_between(m.k_kms,upper,lower,color=[0.8,0.8,0.8],alpha=0.5,label='chi2 weighting')
    plt.title('z={}: alpha={:2.2f}, beta={:2.2f}, sigma_G={:2.2f}, n={:2.2f}, k1={:2.4f}, mean_F={:2.3f} ({:2.3f})'.format(m.z_value,m.alpha,m.beta,m.sigma_G,m.n,m.k1,m.mean_F,tuning.get_mean_F_model(m.z_value)),fontsize=12)
    #plt.axvline(x=0.0152,color='k')
    plt.semilogy()
    plt.semilogx()
    ylim_lower = min(model_Pk_kms) * 0.8
    ylim_upper = max(model_Pk_kms) * 1.2
    plt.ylim(ylim_lower,ylim_upper)
    plt.grid()
    plt.legend(fontsize=12)
    plt.ylabel(r'$P_{1D}$',fontsize=12)
    plt.xlabel(r'$k\ /\ (kms^{-1})^{-1}$',fontsize=12)
    #plt.savefig('Pk1D_z{}.pdf'.format(m.z_value))
    plt.show()

"""
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
"""

"""
#Plot some different parameter values
n_testers = [1.0,1.5,2.03,2.5,3.0]
testers = []
for pixel in pixels:
    for n_tester in n_testers:
        testers += [measure_pixel_segment(pixel,z_value,alpha,beta,sigma_G,n_tester,k1,A0)]
testers = tuning.measurement_set(measurements=testers)
testers = testers.combine_pixels()

plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
k_kms = testers.measurements[0].k_kms
model_Pk_kms = tuning.P1D_z_kms_PD2013(z_value,k_kms)
plt.plot(k_kms,model_Pk_kms,label='DR9 fitting function')
plt.fill_between(k_kms,model_Pk_kms*1.1,model_Pk_kms*0.9,color=[0.5,0.5,0.5],alpha=0.5,label='DR9 +/- 10%')

for m in testers.measurements:
    plt.plot(m.k_kms,m.Pk_kms,label='n = {}'.format(m.n))
plt.title('z={}: alpha={:2.2f}, beta={:2.2f}, sigma_G={:2.2f}, n varied, k1={:2.4f}, mean_F={:2.3f} ({:2.3f})'.format(m.z_value,m.alpha,m.beta,m.sigma_G,m.k1,m.mean_F,tuning.get_mean_F_model(m.z_value)),fontsize=12)
plt.axvline(x=k1,color=[0.25,0.25,0.25])
plt.semilogy()
plt.semilogx()
ylim_lower = min(model_Pk_kms) * 0.8
ylim_upper = max(model_Pk_kms) * 1.2
plt.ylim(ylim_lower,ylim_upper)
plt.grid()
plt.legend(fontsize=12)
plt.ylabel(r'$P_{1D}$',fontsize=12)
plt.xlabel(r'$k\ /\ (kms^{-1})^{-1}$',fontsize=12)
plt.savefig('Pk1D_z{}_ncompare_k1{}.pdf'.format(m.z_value,m.k1))
plt.show()
"""
