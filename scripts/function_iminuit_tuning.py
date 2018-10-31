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

N_files = 32
N_processes = N_files
lambda_min = 3550.0
min_cat_z = 1.8
IVAR_cutoff = 1150.0

#Get the starting values of alpha, beta and sigma_G from file
#Decide which z values we are going to sample at
z_values = [2.0,2.2,2.4,2.6,2.8,3.0,3.2]
z_width = 0.2
colours = ['C0','C1','C2','C3','C4','C5','C6']

cell_size = 0.25 #Mpc/h

max_k = 0.005 #skm-1

#Open up the Gaussian colore files
#base_file_location = '/Users/jfarr/Projects/test_data/test/'
base_file_location = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/process_output_G_hZsmooth_4096_32_sr2.0_bm1_biasG18_picos_newNz_mpz0_seed1003_123_nside16/'
#base_file_location = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v3/v3.0/'
N_side = 16

new_file_structure = '{}/{}/'               #pixel number//100, pixel number
new_filename_structure = '{}-{}-{}.fits'    #file type, nside, pixel number

input_format = 'gaussian_colore'

#get pixels from those directories created by make_master.py
dirs = glob.glob(base_file_location+new_file_structure.format('*','*'))
pixels = []
dirs = dirs[:N_files]
for dir in dirs:
    ending = dir[len(dir)-dir[-2::-1].find('/')-1:-1]
    if ending != 'logs':
        pixels += [int(ending)]
print(pixels)

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


# TODO: want to move this to tuning.py eventually
def measure_pixel_segment(pixel,C0,C1,C2,D0,D1,D2,n,k1,RSD_weights,return_RSD_weights=False):
    start = time.time()
    print('start pixel {} at {}'.format(pixel,time.ctime()))
    location = base_file_location + '/' + new_file_structure.format(pixel//100,pixel)

    # TODO: this is a problem
    measured_SIGMA_G = 1.17

    #We work from the gaussian colore files made in 'pixelise gaussian skewers'.
    gaussian_filename = new_filename_structure.format('gaussian-colore',N_side,pixel)

    #Make a pixel object from it.
    data = simulation_data.SimulationData.get_gaussian_skewers_object(location+gaussian_filename,None,input_format,SIGMA_G=measured_SIGMA_G,IVAR_cutoff=IVAR_cutoff)
    #print('{:3.2f} checkpoint sim_dat'.format(time.time()-start))

    #trim skewers
    data.trim_skewers(lambda_min,min_cat_z,extra_cells=1)

    #extra = 0.1
    #lambda_min_val = lya*(1 + 3.0 - z_width*(1+extra)/2)
    #lambda_max_val = lya*(1 + 3.0 + z_width*(1+extra)/2)
    #data.trim_skewers(lambda_min_val,min_cat_z,lambda_max=lambda_max_val,whole_lambda_range=True)

    #print('Z range',min(data.Z),max(data.Z))
    #print('N_qso',data.N_qso)

    #Expand the C and D parameters into functions of z.
    alpha = C0 * (data.Z**C1) + C2
    beta = np.ones_like(alpha) * 1.65
    sigma_G_required = D0 * (data.Z**D1) + D2

    #Determine the sigma_G to add
    extra_sigma_G = np.sqrt(sigma_G_required**2 - measured_SIGMA_G**2)

    #add small scale fluctuations
    seed = int(str(N_side) + str(pixel))
    #print('seed',seed)
    #print('(n,k1)',n,k1)
    generator = np.random.RandomState(seed)
    data.add_small_scale_gaussian_fluctuations(cell_size,data.Z,extra_sigma_G,generator,amplitude=1.0,white_noise=False,lambda_min=0.0,IVAR_cutoff=lya,n=n,k1=k1,A0=58.6)

    #Copmute the physical skewers
    data.compute_physical_skewers()
    
    #print('Gaussian P1D after SSP seems to match')
    #indices = (data.Z > 2.9)*(data.Z < 3.1)
    #QSOs = np.sum(data.IVAR_rows[:,indices],axis=1)==np.sum(indices)
    #k,pk,_=Pk1D.get_Pk1D(data.GAUSSIAN_DELTA_rows,data.IVAR_rows,data.R,data.Z,z_value=3.0,z_width=0.2,units='km/s')
    #print(k,pk)
    #print('Physical P1D')
    #k,pk,_=Pk1D.get_Pk1D(data.DENSITY_DELTA_rows[:,indices][QSOs,:],data.IVAR_rows[:,indices][QSOs,:],data.R[indices],data.Z[indices],z_value=3.0,z_width=0.2,units='km/s')
    #print(k,pk)
    #print(' ')
    #print('first DD skewer',data.DENSITY_DELTA_rows[0,:10])
    #Update the alpha and beta lengths.
    alpha = C0 * (data.Z**C1) + C2
    beta = np.ones_like(alpha) * 1.65
    sigma_G_required = D0 * (data.Z**D1) + D2
    #print('{:3.2f} checkpoint ssgf'.format(time.time()-start))
    #Compute the tau skewers and add RSDs
    data.compute_tau_skewers(data.lya_absorber,alpha=alpha,beta=beta)

    if return_RSD_weights:
        RSD_weights = data.get_RSD_weights(data.lya_absorber,alpha,beta,thermal=False)
        #print('{:3.2f} checkpoint RSDs measured'.format(time.time()-start))
        return (pixel,RSD_weights)
    else:
        #print('first tau skewer',data.lya_absorber.tau[0,:10])
        RSD_weights = RSD_weights_dict[pixel]
        data.add_RSDs(data.lya_absorber,alpha,beta,thermal=False,weights=RSD_weights)
        #lambda_min_val = lya*(1 + 3.0 - z_width/2)
        #lambda_max_val = lya*(1 + 3.0 + z_width/2)
        #data.trim_skewers(lambda_min_val,min_cat_z,lambda_max=lambda_max_val,whole_lambda_range=True)
        #print('{:3.2f} checkpoint RSDs applied'.format(time.time()-start))
        measurements = []
        for z_value in z_values:
            #print('z={}: alpha={}, sigma_G={}'.format(z_value,np.interp(z_value,data.Z,alpha),np.interp(z_value,data.Z,sigma_G_required)))
            #print('data z:',data.Z)
            #print(pixel,z_value)
            ID = n
            measurement = tuning.function_measurement(ID,z_value,z_width,data.N_qso,n,k1,C0,C1,C2,beta,D0,D1,D2,pixels=[pixel])
            #print('measurement object made')
            measurement.add_mean_F_measurement(data)
            #print('measured mean')
            #print('measuring Pk1D')
            measurement.add_Pk1D_measurement(data)
            #print('measured P1D')
            measurement.add_sigma_F_measurement(data)
            #print('measured sigma')
            measurement.add_mean_F_chi2(eps=0.05)
            measurement.add_Pk1D_chi2(max_k=max_k)
            measurement.add_sigma_F_chi2(eps=0.1)
            measurement.add_total_chi2()
            measurements += [measurement]
            #print(' ')
        #print('{:3.2f} checkpoint measurements'.format(time.time()-start))
        return measurements

#Get the weights for the RSDs
print('getting the weights for the RSDs')
tasks = [(pixel,104.5,-4.62,1.654,54.6,-0.068,-43.81,1.52,0.0166,None,True) for pixel in pixels]

if __name__ == '__main__':
    pool = Pool(processes = N_processes)
    start_time = time.time()
    results = []

    for task in tasks:
        pool.apply_async(measure_pixel_segment,task,callback=log_result,error_callback=log_error)
        #results += [result]

    pool.close()
    pool.join()

RSD_weights_dict = {}
for result in results:
    RSD_weights_dict[result[0]] = result[1]
print('done!')

def f(C0,C1,C2,D0,D1,D2,n,k1,return_measurements=False):

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

    print('looking at params: C=({:2.4f},{:2.4f},{:2.4f}), D=({:2.4f},{:2.4f},{:2.4f}), n={:2.4f}, k1={:2.6f}'.format(C0,C1,C2,D0,D1,D2,n,k1))

    tasks = [(pixel,C0,C1,C2,D0,D1,D2,n,k1,None) for pixel in pixels]
    #print(RSD_weights_dict)
    #manager = multiprocessing.Manager()
    #shared_RSD_weights_dict = manager.dict(RSD_weights_dict)

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

    measurements = []
    for result in results:
        measurements += result

    measurement_set = tuning.measurement_set(measurements=measurements)
    combined_pixels_set = measurement_set.combine_pixels()

    Pk_kms_chi2 = 0.
    mean_F_chi2 = 0.
    overall_chi2 = 0.

    sigma_F_chi2 = 0.

    print('checkpoint')
    for m in combined_pixels_set.measurements:

        m.add_mean_F_chi2(eps=0.05)
        m.add_Pk1D_chi2(max_k=max_k)#,denom="npower")
        m.add_sigma_F_chi2(eps=0.05)
        m.add_total_chi2()
        Pk_kms_chi2 += m.Pk_kms_chi2
        mean_F_chi2 += m.mean_F_chi2
        overall_chi2 += m.total_chi2
        #print(m.z_value)
        #print(m.mean_F)
        #print(m.k_kms)
        #print(m.Pk_kms)

        #print(m.get_details())
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

#Parameters defined by:
#   alpha = C0 * (Z**C1) + C2
#   beta = np.ones_like(alpha) * 1.65
#   sigma_G_required = D0 * (Z**D1) + D2

fix_all = False

a_kwargs = {'C0' : 104.5,     'error_C0' : 1.0,  'fix_C0' : fix_all|False, #'limit_C0' : (0., 20.),
            'C1' : -4.632,     'error_C1' : 0.05,  'fix_C1' : fix_all|False, #'limit_C1' : (0., 20.),
            'C2' : 1.654,     'error_C2' : 0.05,  'fix_C2' : fix_all|False, #'limit_C2' : (0., 20.),
            }

sG_kwargs = {'D0' : 54.6,     'error_D0' : 1.0,  'fix_D0' : fix_all|False, #'limit_D0' : (0., 20.),
             'D1' : -0.068,     'error_D1' : 0.01,  'fix_D1' : fix_all|False, #'limit_D1' : (0., 20.),
             'D2' : -43.81,     'error_D2' : 0.5,  'fix_D2' : fix_all|False, #'limit_D2' : (0., 20.),
             }

s_kwargs = {'n'  : 1.52,     'error_n' : 0.05,   'limit_n' : (0., 10.),   'fix_n' : fix_all|False,
            'k1' : 0.0166,   'error_k1' : 0.0005,'limit_k1' : (0., 0.1),  'fix_k1' : fix_all|False,
            }

other_kwargs = {'return_measurements'  : False,    'fix_return_measurements' : True,
            }

minuit = Minuit(f,**a_kwargs,**sG_kwargs,**s_kwargs,**other_kwargs)

minuit.print_param()
minuit.migrad() #ncall=20
minuit.print_param()

C0 = minuit.values['C0']
C1 = minuit.values['C1']
C2 = minuit.values['C2']
#beta = minuit.values['beta']
D0 = minuit.values['D0']
D1 = minuit.values['D1']
D2 = minuit.values['D2']
n = minuit.values['n']
k1 = minuit.values['k1']

print(minuit.values)

#Want to do a final run here
#final_measurements = []
#for pixel in pixels:
#    final_measurements += [measure_pixel_segment(pixel,C0,C1,C2,D0,D1,D2,n,k1,RSD_weights_dict[pixel],return_RSD_weights=False)]
#final_measurements = tuning.measurement_set(measurements=final_measurements)
#final_measurements = final_measurements.combine_pixels()

final_chi2,final_measurements = f(C0,C1,C2,D0,D1,D2,n,k1,return_measurements=True)


#Plot a graph of mean F with redshift
def plot_mean_F_values(m_set,wait=False):
    mean_F = []
    for m in m_set.measurements:
        mean_F += [m.mean_F]
    plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
    plt.plot(z_values,mean_F,label='mean flux',marker='o')
    model_mean_F = tuning.get_mean_F_model(np.array(z_values))
    plt.plot(z_values,model_mean_F,label='model',marker='o')
    plt.fill_between(z_values,model_mean_F*1.1,model_mean_F*0.9,color=[0.5,0.5,0.5],alpha=0.5)
    plt.grid()
    plt.legend()
    if wait:
        plt.show()
    else:
        plt.pause(0.01)
    return

#Plot the power spectra
def plot_P1D_values(m_set,wait=False):
    plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
    for m in m_set.measurements:
        colour = colours[np.searchsorted(z_values,m.z_value)]
        plt.plot(m.k_kms,m.Pk_kms,label='z={}'.format(m.z_value),c=colour)
        model_Pk_kms = tuning.P1D_z_kms_PD2013(m.z_value,m.k_kms)
        plt.plot(m.k_kms,model_Pk_kms,label='DR9 z={}'.format(m.z_value),c=colour,linestyle=':')
        m.add_Pk1D_chi2(max_k=max_k)#,denom="npower")
        eps = m.Pk_kms_chi2_eps
        plt.fill_between(m.k_kms,model_Pk_kms*1.1,model_Pk_kms*0.9,color=[0.5,0.5,0.5],alpha=0.5)#,label='DR9 +/- 10%')
        lower = np.maximum(np.ones_like(model_Pk_kms)*10**(-6),model_Pk_kms * (1. - eps))
        upper = model_Pk_kms * (1. + eps)
        #plt.plot(m.k_kms,upper,c='k',linestyle='dashed')
        #plt.plot(m.k_kms,lower,c='k',linestyle='dashed')
        #plt.fill_between(m.k_kms,upper,lower,color=[0.8,0.8,0.8],alpha=0.5,label='chi2 weighting')
    plt.title('z={}: C=({:2.4f},{:2.4f},{:2.4f}), D=({:2.4f},{:2.4f},{:2.4f}), n={:2.4f}, k1={:2.6f}'.format(m.z_value,m.C0,m.C1,m.C2,m.D0,m.D1,m.D2,m.n,m.k1),fontsize=12)
    plt.axvline(x=max_k,color='k')
    plt.semilogy()
    plt.semilogx()
    ylim_lower = min(model_Pk_kms) * 0.8
    ylim_upper = max(model_Pk_kms) * 1.2
    plt.ylim(ylim_lower,ylim_upper)
    plt.xlim(xmax=0.02)
    plt.grid()
    plt.legend(fontsize=12)
    plt.ylabel(r'$P_{1D}$',fontsize=12)
    plt.xlabel(r'$k\ /\ (kms^{-1})^{-1}$',fontsize=12)
    #plt.savefig('Pk1D_z{}.pdf'.format(m.z_value))
    if wait:
        plt.show()
    else:
        plt.pause(0.01)
    return

plot_mean_F_values(final_measurements,wait=True)
plot_P1D_values(final_measurements,wait=True)

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
