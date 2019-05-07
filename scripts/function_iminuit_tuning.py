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

lya = utils.lya_rest

N_files = 32
N_processes = np.min((64,N_files))
lambda_min = 3550.0
min_cat_z = 1.8
IVAR_cutoff = 1150.0
global_seed = 123
add_ssf = True
lambda_buffer = 100. #Angstroms

#Choose parameter values.
eps_Pk1D = 0.1
eps_mean_F = 0.025
eps_bias_delta = 0.025
eps_bias_eta = 10**6#0.025
d_delta = 10.**-3
d_eta = 10**-9

#Choose tuning parameter initial values.
initial_C0 = 1.4665375029450034
initial_C1 = 4.5
initial_C2 = 0.0
initial_beta = 1.65
initial_D0 = 6.103687317112405
initial_D1 = 0.32156149272713025
initial_D2 = 0.0
initial_n = 0.8322104177553062
initial_k1 = 0.017492626170643323
initial_R = 25.0 #kms-1
initial_vb = 1.0

#Choose parameters to fix.
fix_all = False
fix_C0 = False
fix_C1 = True
fix_C2 = True
fix_beta = True
fix_D0 = False
fix_D1 = False
fix_D2 = True
fix_n = False
fix_k1 = False
fix_R = True
fix_vb = True

#Admin options
k_plot_max = 0.1
show_plots = False
save_plots = True
suffix = '_with_bias_{}RSD_vel{}_a{}_b{}'.format('new','NGP','free',initial_beta)
save_tuning = True
overwrite_tuning = True
tuning_filename = 'input_files/tuning_data' + suffix + '.fits'

#Get the starting values of alpha, beta and sigma_G from file
#Decide which z values we are going to sample at
z_values = [2.0,2.2,2.4,2.6,2.8,3.0,3.2]
z_width = 0.1
colours = ['C0','C1','C2','C3','C4','C5','C6']
cell_size = 0.25 #Mpc/h
max_k = 0.01 #skm-1

#Open up the Gaussian colore files
base_file_location = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v5/v5.0.0/'
N_side = 16

new_file_structure = '{}/{}/'               #pixel number//100, pixel number
new_filename_structure = '{}-{}-{}.fits'    #file type, nside, pixel number

input_format = 'gaussian_colore'

#Function to produce a parameter as a function of z.
def get_parameter(z,A0,A1,A2,z0=3.0):
    x = (1 + z)/(1 + z0)
    alpha = np.exp(utils.quadratic_log(x,A0,A1,A2))
    return alpha

#get pixels from those directories created by make_master.py
dirs = glob.glob(base_file_location+new_file_structure.format('*','*'))
pixels = []
dirs = dirs[:N_files]
for dir in dirs:
    ending = dir[len(dir)-dir[-2::-1].find('/')-1:-1]
    if ending != 'logs':
        pixels += [int(ending)]
pixels = list(range(N_files))

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
def measure_pixel_segment(pixel,C0,C1,C2,beta_value,D0,D1,D2,n,k1,R_kms,vel_boost,RSD_weights,prep=False):
    t = time.time()
    seed = int(pixel * 10**5 + global_seed)

    #print('start pixel {} at {}'.format(pixel,time.ctime()))
    location = base_file_location + '/' + new_file_structure.format(pixel//100,pixel)

    #We work from the gaussian colore files made in 'pixelise gaussian skewers'.
    gaussian_filename = new_filename_structure.format('gaussian-colore',N_side,pixel)

    #Make a pixel object from it.
    data = simulation_data.SimulationData.get_gaussian_skewers_object(location+gaussian_filename,None,input_format,IVAR_cutoff=IVAR_cutoff)
    #print('{:3.2f} checkpoint sim_dat'.format(time.time()-t))
    t = time.time()

    #Scale the RSD skewers.
    data.VEL_rows *= vel_boost

    transformation = tuning.transformation()
    def f_tau0_z(z):
        return get_parameter(z,C0,C1,C2)
    def f_texp_z(z):
        return get_parameter(z,beta_value,0.,0.)
    def f_seps_z(z):
        return get_parameter(z,D0,D1,D2)
    transformation.add_parameters_from_functions(f_tau0_z,f_texp_z,f_seps_z)
    data.transformation = transformation

    #trim skewers to the minimal length
    extra = 0.1
    z_lower_cut = np.min(z_values) - z_width*(1+extra)/2.
    z_upper_cut = np.max(z_values) + z_width*(1+extra)/2.
    lambda_min_val = np.min([lambda_min,lya*(1 + z_lower_cut)])
    lambda_max_val = lya*(1 + z_upper_cut)
    data.trim_skewers(lambda_min_val,min_cat_z,lambda_max=lambda_max_val,whole_lambda_range=False)

    #lambda_min_val = lambda_min-100.
    #data.trim_skewers(lambda_min_val,min_cat_z,extra_cells=1)

    #add small scale fluctuations
    if add_ssf:
        generator = np.random.RandomState(seed)
        data.add_small_scale_gaussian_fluctuations(cell_size,generator,white_noise=False,lambda_min=0.0,IVAR_cutoff=IVAR_cutoff,n=n,k1=k1,R_kms=R_kms)

        #print('{:3.2f} checkpoint extra power'.format(time.time()-t))
        t = time.time()

    #Copmute the physical skewers
    data.compute_physical_skewers()

    #Way to get exact same alpha as in make_transmission
    #z = np.linspace(1.6,4.0,2401)
    #alpha = get_alpha(z,C0,C1,C2)
    #beta = np.ones_like(alpha) * beta_value
    #extra_sigma_G = get_sigma_G(z,D0,D1,D2)
    #alpha = np.exp(np.interp(np.log(data.Z),np.log(z),np.log(alpha)))
    #beta = np.exp(np.interp(np.log(data.Z),np.log(z),np.log(beta)))
    #extra_sigma_G = np.exp(np.interp(np.log(data.Z),np.log(z),np.log(extra_sigma_G)))

    """
    #HACK
    manual_cells = data.Z<2.5
    transfer_cell = np.searchsorted(data.Z,2.5)
    grad = (alpha[transfer_cell] - alpha[transfer_cell + 1])/(data.Z[transfer_cell + 1] - data.Z[transfer_cell])
    manual_alphas = alpha[transfer_cell] + grad*(data.Z[transfer_cell] - data.Z)
    alpha = alpha - manual_cells * (alpha - manual_alphas)
    """

    #Compute the tau skewers and add RSDs
    data.compute_tau_skewers(data.lya_absorber)
    #print('{:3.2f} checkpoint tau'.format(time.time()-t))
    t = time.time()

    if prep:
        RSD_weights = data.get_RSD_weights(thermal=False)
        #print(pixel,'{:3.2f} checkpoint RSD weights measured'.format(time.time()-t))
        t = time.time()

        b_eta_weights_dict = data.get_bias_eta_RSD_weights(z_values,d=d_eta,z_width=z_width,lambda_buffer=lambda_buffer)

        #print(pixel,'{:3.2f} checkpoint b_eta weights measured'.format(time.time()-t))
        t = time.time()

        return (pixel,RSD_weights,b_eta_weights_dict)
    else:
        RSD_weights = RSD_weights_dict[pixel]
        bias_eta_weights = bias_eta_weights_dict[pixel]
        data.add_all_RSDs(thermal=False,weights=RSD_weights)

        #print('{:3.2f} checkpoint RSDs'.format(time.time()-t))
        t = time.time()

        measurements = []
        times_m = np.zeros(6)
        for z_value in z_values:
            ID = n
            t_m = time.time()
            measurement = tuning.function_measurement(ID,z_value,z_width,data.N_qso,n,k1,C0,C1,C2,beta_value,D0,D1,D2,pixels=[pixel])
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
            measurement.add_bias_eta_measurement(data,d=d_eta,weights_dict=bias_eta_weights,lambda_buffer=lambda_buffer)
            times_m[5] += time.time() - t_m
            measurements += [measurement]

        #print('{:3.2f} checkpoint measurements'.format(time.time()-t))
        #print('--> measurement_times: {:3.2f}, {:3.2f}, {:3.2f}, {:3.2f}, {:3.2f}, {:3.2f}'.format(times_m[0],times_m[1],times_m[2],times_m[3],times_m[4],times_m[5]))

        return measurements

#Pre-prep for future processings by getting RSD maps and independent skewers
print('producing preparatory data (RSD maps)')
#tasks = [(pixel,104.5,-4.62,1.654,54.6,-0.068,-43.81,1.52,0.0166,None,True) for pixel in pixels]
tasks = [(pixel,initial_C0,initial_C1,initial_C2,initial_beta,initial_D0,initial_D1,initial_D2,initial_n,initial_k1,initial_R,initial_vb,None,True) for pixel in pixels]

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
bias_eta_weights_dict = {}
for result in results:
    RSD_weights_dict[result[0]] = result[1]
    bias_eta_weights_dict[result[0]] = result[2]
print('done!')

def f(C0,C1,C2,beta,D0,D1,D2,n,k1,R,vb,return_measurements=False):

    ################################################################################
    """
    #Define the multiprocessing tracking functions
    """

    #Define a progress-tracking function.
    def log_result(retval):
        results.append(retval)
        return retval

    #Define an error-tracking function.
    def log_error(retval):
        print('Error:',retval)

    ################################################################################

    print('starting at',time.ctime())
    print('looking at params: C=({:2.4f},{:2.4f},{:2.4f}), beta={:1.2f},  D=({:2.4f},{:2.4f},{:2.4f}), n={:2.4f}, k1={:2.6f}'.format(C0,C1,C2,beta,D0,D1,D2,n,k1))

    tasks = [(pixel,C0,C1,C2,beta,D0,D1,D2,n,k1,R,vb,None) for pixel in pixels]

    #Run the multiprocessing pool
    if __name__ == '__main__':
        pool = Pool(processes = N_processes)
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
        m.add_bias_eta_chi2(eps=eps_bias_eta)

        Pk_kms_chi2 += m.Pk_kms_chi2
        mean_F_chi2 += m.mean_F_chi2
        bias_delta_chi2 += m.bias_delta_chi2
        bias_eta_chi2 += m.bias_eta_chi2
        #print('z =',m.z_value,'number of k values =',m.k_kms.shape,'number of k values < max k =',np.sum(m.k_kms<max_k))

    chi2 = mean_F_chi2 + Pk_kms_chi2 + bias_delta_chi2 + bias_eta_chi2
    log_text = 'chi2: Pk {:2.4f}, mean F {:2.4f}, b_del {:2.4f}, b_eta {:2.4f}, overall {:2.4f}'.format(Pk_kms_chi2,mean_F_chi2,bias_delta_chi2,bias_eta_chi2,chi2)

    print(log_text)
    print(' ')

    with open("parameter_log.txt","a") as f:
        txt = str(time.ctime()+'\n')
        txt += 'C0:{}, C1:{}, C2:{}, D0:{}, D1:{}, D2:{}, n:{}, k1:{}, beta:{}\n'.format(C0,C1,C2,D0,D1,D2,n,k1,beta)
        txt += log_text + '\n\n'
        f.write(txt)
        f.close()
        best = chi2

    if return_measurements:
        return chi2, combined_pixels_set
    else:
        return chi2

#Parameters defined by:
#   alpha = C0 * (Z**C1) + C2
#   beta = np.ones_like(alpha) * beta_value
#   sigma_G_required = D0 * (Z**D1) + D2


a_kwargs = {'C0' : initial_C0,      'error_C0' : 1.0,   'fix_C0' : fix_all|fix_C0,     'limit_C0' : (0., 100.),
            'C1' : initial_C1,      'error_C1' : 1.0,   'fix_C1' : fix_all|fix_C1,     #'limit_C1' : (0., 20.),
            'C2' : initial_C2,      'error_C2' : 1.0,   'fix_C2' : fix_all|fix_C2,     #'limit_C2' : (0., 20.),
            }

b_kwargs = {'beta' : initial_beta,  'error_beta' : 1.0, 'fix_beta' : fix_all|fix_beta, 'limit_beta' : (0.,5.)
            }

sG_kwargs = {'D0' : initial_D0,     'error_D0' : 1.0,   'fix_D0' : fix_all|fix_D0,     'limit_D0' : (0., 100.),
             'D1' : initial_D1,     'error_D1' : 0.2,   'fix_D1' : fix_all|fix_D1,     #'limit_D1' : (0., 20.),
             'D2' : initial_D2,     'error_D2' : 1.0,   'fix_D2' : fix_all|fix_D2,     #'limit_D2' : (0., 20.),
             }

s_kwargs = {'n'  : initial_n,       'error_n' : 1.0,    'fix_n' : fix_all|fix_n,       'limit_n' : (-2., 10.),
            'k1' : initial_k1,      'error_k1' : 0.001, 'fix_k1' : fix_all|fix_k1,     'limit_k1' : (0., 0.1),
            }

other_kwargs = {'R'  : initial_R,    'error_R' : 1.0,   'fix_R' : fix_all|fix_R,       'limit_R' : (0., 1000.),
                'vb' : initial_vb,   'error_vb' : 0.1,  'fix_vb' : fix_all|fix_vb,     'limit_vb' : (0., 2.0),
                'return_measurements'  : False,    'fix_return_measurements' : True,
                'errordef'             : 1,
                }

minuit = Minuit(f,**a_kwargs,**b_kwargs,**sG_kwargs,**s_kwargs,**other_kwargs)

if not fix_all:
    minuit.print_param()
    minuit.migrad() #ncall=20
    minuit.print_param()

C0 = minuit.values['C0']
C1 = minuit.values['C1']
C2 = minuit.values['C2']
beta = minuit.values['beta']
D0 = minuit.values['D0']
D1 = minuit.values['D1']
D2 = minuit.values['D2']
n = minuit.values['n']
k1 = minuit.values['k1']
R = minuit.values['R']
vb = minuit.values['vb']

print(minuit.values)

def save_tuning_file(filename,overwrite=False):
    z = np.linspace(1.6,4.0,2401)
    alpha_arr = get_parameter(z,C0,C1,C2)
    beta_arr = np.ones(alpha_arr.shape)*beta
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

    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)

    dtype = [('z', 'f4'), ('alpha', 'f4'), ('beta', 'f4'), ('sigma_G', 'f4')]
    data = np.array(list(zip(z,alpha_arr,beta_arr,sigma_G_arr)),dtype=dtype)
    hdu_tuning = fits.BinTableHDU(data,header=header,name='TUNING')

    hdulist = fits.HDUList([prihdu,hdu_tuning])
    hdulist.writeto(filename,overwrite=overwrite)
    hdulist.close

    return

if save_tuning:
    save_tuning_file(tuning_filename,overwrite=overwrite_tuning)

#Do a final run to get the measurements.
final_chi2,final_measurements = f(C0,C1,C2,beta,D0,D1,D2,n,k1,R,vb,return_measurements=True)

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
        model = tuning.get_mean_F_model(np.array(z_values))
    elif data_type == 'bias_delta':
        text = r'$b_\delta$'
        model = tuning.get_bias_delta_model(np.array(z_values))
    elif data_type == 'bias_eta':
        text = r'$b_\eta$'
        model = tuning.get_bias_eta_model(np.array(z_values))

    plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
    plt.plot(data['z'],data['data'],label=text+' mocks',marker='o')

    plt.plot(z_values,model,label=text+' model',marker='x',c=[0.5,0.5,0.5])
    plt.fill_between(z_values,model*1.1,model*0.9,color=[0.5,0.5,0.5],alpha=0.3)
    plt.grid()
    plt.legend()

    z = data['z']
    min_z = np.min(z) - (np.max(z) - np.min(z))*0.1
    max_z = np.max(z) + (np.max(z) - np.min(z))*0.1
    plt.xlim(min_z,max_z)
    plt.xlabel(r'$z$',fontsize=12)
    plt.ylabel(text,fontsize=12)
    if save_plot:
        plt.savefig(data_type + suffix + '.pdf')
    if show_plot:
        plt.show()
    return

#Plot the power spectra
def plot_P1D_values(m_set,show_plot=True,save_plot=False):
    plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
    max_power_plot = 0.
    min_power_plot = 10.**10
    for m in m_set.measurements:
        colour = colours[np.searchsorted(z_values,m.z_value)]
        plt.plot(m.k_kms,m.Pk_kms,label=r'$z={}$'.format(m.z_value),c=colour)
        model_Pk_kms = tuning.P1D_z_kms_PD2013(m.z_value,m.k_kms)
        plt.plot(m.k_kms,model_Pk_kms,c=colour,linestyle=':')#,label='z={} DR9'.format(m.z_value))
        m.add_Pk1D_chi2(max_k=max_k,denom="npower_cutoff")
        eps = m.Pk_kms_chi2_eps
        #plt.plot(m.k_kms,model_Pk_kms*0.9,color=[0.5,0.5,0.5],alpha=0.5)
        #plt.plot(m.k_kms,model_Pk_kms*1.1,color=[0.5,0.5,0.5],alpha=0.5)
        plt.fill_between(m.k_kms,model_Pk_kms*1.1,model_Pk_kms*0.9,color=[0.5,0.5,0.5],alpha=0.3)#,label='DR9 +/- 10%')
        lower = np.maximum(np.ones_like(model_Pk_kms)*10**(-6),model_Pk_kms * (1. - eps))
        upper = model_Pk_kms * (1. + eps)
        plt.plot(m.k_kms,upper,c='k',linestyle='dashed')
        plt.plot(m.k_kms,lower,c='k',linestyle='dashed')
        #plt.fill_between(m.k_kms,upper,lower,color=[0.8,0.8,0.8],alpha=0.3)
        max_power_plot = np.max((max_power_plot,np.max(model_Pk_kms[m.k_kms<k_plot_max])))
        min_power_plot = np.min((min_power_plot,np.min(model_Pk_kms[m.k_kms<k_plot_max])))
    plt.title('C=({:2.4f},{:2.4f},{:2.4f}), D=({:2.4f},{:2.4f},{:2.4f}), n={:2.4f}, k1={:2.6f}'.format(m.C0,m.C1,m.C2,m.D0,m.D1,m.D2,m.n,m.k1),fontsize=12)
    plt.axvline(x=max_k,color='k')
    plt.semilogy()
    plt.semilogx()
    ylim_lower = min_power_plot * 0.8
    ylim_upper = max_power_plot * 1.2
    plt.ylim(ylim_lower,ylim_upper)
    plt.xlim(xmax=k_plot_max)
    plt.grid()
    plt.legend(fontsize=12)
    plt.ylabel(r'$P_{1D}$',fontsize=12)
    plt.xlabel(r'$k\ /\ (kms^{-1})^{-1}$',fontsize=12)
    if save_plot:
        plt.savefig('Pk1D' + suffix + '.pdf')
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
        plt.savefig('parameters' + suffix + '.pdf')
    if show_plot:
        plt.show()
    return

plot_scalar_measurement_values(final_measurements,data_type='mean_F',show_plot=show_plots,save_plot=save_plots)
plot_scalar_measurement_values(final_measurements,data_type='bias_delta',show_plot=show_plots,save_plot=save_plots)
plot_scalar_measurement_values(final_measurements,data_type='bias_eta',show_plot=show_plots,save_plot=save_plots)
plot_P1D_values(final_measurements,show_plot=show_plots,save_plot=save_plots)
plot_parameter_values(minuit,z_min=min(z_values)-z_width/2.,z_max=max(z_values)+z_width/2.,show_plot=show_plots,save_plot=save_plots)
