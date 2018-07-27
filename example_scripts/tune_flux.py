import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
import subprocess
from multiprocessing import Pool
import multiprocessing
import time
import glob

import pixelise
import Pk1D
import tuning
import independent
import general

lya = 1215.67

N_processes = 4
lambda_min = 3550.0
min_cat_z = 1.8
IVAR_cutoff = 1150.0

#Get the starting values of alpha, beta and sigma_G from file
#Decide which z values we are going to tune
z_values = [2.5]
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
for dir in dirs[:1]:
    pixels += [int(dir[len(dir)-dir[-2::-1].find('/')-1:-1])]
#pixels=[0]

################################################################################

"""
Define the multiprocessing tracking functions
"""

#Define a progress-tracking function.
def log_result(retval):

    results.append(retval)
    N_complete = len(results)
    N_tasks = len(tasks)

    general.progress_bar(N_complete,N_tasks,start_time)

#Define an error-tracking function.
def log_error(retval):
    print('Error:',retval)

################################################################################

#load initial tuning data
tuning_file = 'input_files/tune_small_scale_fluctuations.fits'
h = fits.open(tuning_file)
tune_small_scale_fluctuations = h['DATA'].data
h.close()
tuning_z_values = tune_small_scale_fluctuations['z']
desired_sigma_G_values = tune_small_scale_fluctuations['sigma_G']
desired_mean_F = tune_small_scale_fluctuations['mean_F']
alphas = tune_small_scale_fluctuations['alpha']

beta_default = 1.65
betas = beta_default*np.ones(tuning_z_values.shape[0])

n_default = 0.7
ns = n_default*np.ones(tuning_z_values.shape[0])

k1_default = 0.001
k1s = k1_default*np.ones(tuning_z_values.shape[0])

#Get the values of alpha, beta, sigma_G for the z_values
alpha_values = []
beta_values = []
sigma_G_values = []
n_values = []
k1_values = []

for z_value in z_values:
    alpha_values += [np.interp(z_value,tuning_z_values,alphas)]
    beta_values += [np.interp(z_value,tuning_z_values,betas)]
    sigma_G_values += [np.interp(z_value,tuning_z_values,desired_sigma_G_values)]
    n_values += [np.interp(z_value,tuning_z_values,ns)]
    k1_values += [np.interp(z_value,tuning_z_values,k1s)]


################################################################################
"""
#Set the parameter values that we will iterate over
import itertools
for i,z_value in enumerate(z_values):
    a = alpha_values[i] * multipliers
    b = [beta_values[i]]
    sG = sigma_G_values[i] * multipliers

"""
#s_multipliers = np.linspace(0.7,1.3,5)
t_multipliers = np.linspace(0.8,1.2,1)
n_multipliers = np.linspace(0.8,1.2,1)
k1_multipliers = np.linspace(0.8,1.2,1)

#Extract the values of parameters to optimise over
parameters_list = []

lookup = {}

import itertools
for i,z_value in enumerate(z_values):
    a = [alpha_values[i]] * t_multipliers
    b = [beta_values[i]] * t_multipliers
    sG = [sigma_G_values[i]] * t_multipliers

    n = [n_values[i]] * n_multipliers
    k1 = [k1_values[i]] * k1_multipliers

    #parameters_list += list(itertools.product([z_value],a,b,sG,n,k1))

    s_parameter_values_list = list(itertools.product(n,k1))
    t_parameter_values_list = list(itertools.product(a,b,sG))

    # TODO: need a better ID number of some kind
    z_dict = {}
    ID = 0
    for s_parameter_values in s_parameter_values_list:
        for t_parameter_values in t_parameter_values_list:
            s_parameters_dict = {'n':s_parameter_values[0],'k1':s_parameter_values[1]}
            t_parameters_dict = {'alpha':t_parameter_values[0],'beta':t_parameter_values[1],'sigma_G':t_parameter_values[2]}

            Pk_dict = {'k_kms':None,'Pk_kms':None,'var_kms':None}
            results_dict = {'mean_F':None,'Pk1D':Pk_dict,'A_F':None,'B_F':None,'bias':None,'beta':None}
            errors_dict = {'mean_F':None,'Pk1D':None,'A_F':None,'B_F':None,'bias':None,'beta':None}
            ID_dict = {'s_parameters':s_parameters_dict,'t_parameters':t_parameters_dict,'results':results_dict,'errors':errors_dict}
            z_dict = {**z_dict,**{ID:ID_dict}}
            ID += 1


    #data_type = [('alpha','f4'),('beta','f4'),('sigma_G','f4'),('n','f4'),('k1','f4')]
    #parameters = np.array(parameter_values,dtype=data_type)

    lookup = {**lookup,**{z_value:z_dict}}

ID_list = list(range(ID))

def measure_pixel_segment(pixel,z_value,ID,lookup):
    #print('start',z_value)

    n = lookup[z_value][ID]['s_parameters']['n']
    k1 = lookup[z_value][ID]['s_parameters']['k1']

    alpha = lookup[z_value][ID]['t_parameters']['alpha']
    beta = lookup[z_value][ID]['t_parameters']['beta']
    sigma_G_required = lookup[z_value][ID]['t_parameters']['sigma_G']

    location = base_file_location + '/' + new_file_structure.format(pixel//100,pixel)
    mean_F_data = np.array(list(zip(tuning_z_values,desired_mean_F)))

    # TODO: this is a problem
    measured_SIGMA_G = 1.17

    #We work from the gaussian colore files made in 'pixelise gaussian skewers'.
    gaussian_filename = new_filename_structure.format('gaussian-colore',N_side,pixel)

    #Make a pixel object from it.
    data = pixelise.simulation_data.get_gaussian_skewers_object(location+gaussian_filename,None,input_format,SIGMA_G=measured_SIGMA_G,IVAR_cutoff=IVAR_cutoff)

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
    data.add_small_scale_gaussian_fluctuations(cell_size,data.Z,np.ones(data.Z.shape[0])*extra_sigma_G,generator,amplitude=1.0,white_noise=False,lambda_min=0.0,IVAR_cutoff=lya,n=n,k1=k1) #n=0.7, k1=0.001 default

    #Convert to flux
    data.compute_physical_skewers()
    data.compute_tau_skewers(alpha=np.ones(data.Z.shape[0])*alpha,beta=beta)
    data.add_RSDs(np.ones(data.Z.shape[0])*alpha,beta,thermal=False)
    data.compute_flux_skewers()

    #Trim the skewers again to get rid of the additional cells
    lambda_min_val = lya*(1 + z_value - z_width/2)
    lambda_max_val = lya*(1 + z_value + z_width/2)
    data.trim_skewers(lambda_min_val,min_cat_z,lambda_max=lambda_max_val,whole_lambda_range=True)

    #Convert back to small cells for cf measurement
    #need a new function to merge cells back together

    measurement = tuning.measurement(ID,z_value,z_width,data.N_qso,n,k1,alpha,beta,sigma_G_required,pixels=[pixel])
    #print('measure mean_F',z_value)
    measurement.add_mean_F_measurement(data)
    #print('measure Pk1D',z_value)
    measurement.add_Pk1D_measurement(data)
    #print('measure chi2s',z_value)
    measurement.add_mean_F_chi2(eps=0.05)
    measurement.add_Pk1D_chi2(max_k=max_k)
    measurement.add_total_chi2()

    return measurement


tasks = [(pixel,z_value,ID,lookup) for ID in ID_list for pixel in pixels for z_value in z_values]

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = N_processes)
    start_time = time.time()
    results = []

    for task in tasks:
        pool.apply_async(measure_pixel_segment,task,callback=log_result,error_callback=log_error)

    pool.close()
    pool.join()

measurement_set = tuning.measurement_set(measurements=results)

print('originally number measuremts is:',len(measurement_set.measurements))
combined_set = measurement_set.combine_pixels()
print('after combining, number measuremts is:',len(combined_set.measurements))

#Fit model
def get_model_Pk_kms(k_kms,A_F,B_F):
    return tuning.P1D_z_kms_PD2013(k_kms,z_value,A_F=A_F,B_F=B_F)

measurement_set.save('measurements.fits',existing='combine')
optimal_measurements = measurement_set.optimize_s_parameters(plot_optimal=True)


"""

s_parameters_chi2 = []
best_measurements = np.zeros((len(s_parameter_values_list),len(z_values)))
for i,s_parameter_values in enumerate(s_parameter_values_list):
    n = s_parameter_values[0]
    k1 = s_parameter_values[1]
    s_set = measurement_set.s_filter(n,k1)
    total_chi2 = 0
    for j,z_value in enumerate(z_values):
        z_s_set = s_set.z_filter(z_value)
        best = z_s_set.get_best_measurement()
        max_j = np.searchsorted(best.k_kms,max_k)
        fit = curve_fit(get_model_Pk_kms,best.k_kms[:max_j],best.Pk_kms[:max_j],p0=(0.064,3.55))
        best_measurements[i,j] = best
        total_chi2 += best.total_chi2
    s_parameters_chi2 += [total_chi2]

#This is very basic. Improve
i = np.argmin(s_parameters_chi2)
print('best overall s parameters are:',s_parameter_values_list[i])
n_best = s_parameter_values_list[i][0]
k1_best = s_parameter_values_list[i][1]

plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
for j in range(len(z_values)):
    m = best_measurements[i,j]




for z_value in z_values:
    print('z value: {}'.format(z_value))
    z_set = measurement_set.z_filter(z_value)
    for s_parameter_values in s_parameter_values_list:
        n = s_parameter_values[0]
        k1 = s_parameter_values[1]
        z_s_set = z_set.s_filter(n,k1)
        for m in z_s_set.measurements:
            m.add_Pk1D_chi2(max_k=max_k)
            m.add_mean_F_chi2(eps=0.01)
        best = z_s_set.get_best_measurement()
        max_j = np.searchsorted(best.k_kms,max_k)
        fit = curve_fit(get_model_Pk_kms,best.k_kms[:max_j],best.Pk_kms[:max_j],p0=(0.064,3.55))

        print('alpha = {:2.2f}'.format(best.alpha))
        print('beta = {:2.2f}'.format(best.beta))
        print('sigma_G = {:2.2f}'.format(best.sigma_G))
        print('n = {:2.2f}'.format(best.n))
        print('k1 = {:2.4f}'.format(best.k1))
        print(' ')
        print('A_F = {:2.2f}, error is {:3.0%}'.format(fit[0][0],(fit[0][0]/0.064 - 1)))
        print('B_F = {:2.2f}, error is {:3.0%}'.format(fit[0][1],(fit[0][1]/3.55 - 1)))
        print(' ')
        print('mean F = {:2.2f}, chi2 is {:2.2f}'.format(best.mean_F,best.mean_F_chi2))
        print('Pk1D chi2 is {:2.2f}'.format(best.Pk_kms_chi2))

        plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
        plt.title('z_value={:2.2f}: alpha={:2.2f}, beta={:2.2f}, sG={:2.2f}, n={:2.2f}, k1={:2.4f}'.format(z_value,best.alpha,best.beta,best.sigma_G,best.n,best.k1))
        ##Need to update the errors in here
        plt.errorbar(best.k_kms,best.Pk_kms,yerr=np.zeros(best.k_kms.shape),fmt='o',label='measured',color='orange')
        plt.plot(best.k_kms,get_model_Pk_kms(best.k_kms,fit[0][0],fit[0][1]),label='model fit: A_F={:2.2f}, B_F={:2.2f}'.format(fit[0][0],fit[0][1]),color='b')
        plt.plot(best.k_kms,get_model_Pk_kms(best.k_kms,0.064,3.55),label='model default: A_F={:2.2f}, B_F={:2.2f}'.format(0.064,3.55),color=(0.5,0.5,0.5))
        plt.axvline(x=max_k,c=(0.,0.,0.),label='max fitting k value')
        #plt.fill_between(k_kms,0.9*model_Pk_kms_default,1.1*model_Pk_kms_default,label='model default +/- 10%',color=(0.5,0.5,0.5),alpha=0.5)
        #plt.plot(k_kms,independent.power_kms(z_value,k_kms,cell_size*general.get_dkms_dhMpc(z_value),False),label='added')
        plt.semilogy()
        plt.semilogx()
        plt.ylabel('Pk1D')
        plt.xlabel('k / kms-1')
        plt.legend()
        plt.grid()
        #plt.savefig('Pk1D_abs_slope1.5_alpha{:2.2f}_beta{:2.2f}_sG{:2.2f}.pdf'.format(alpha,beta,sigma_G_required))
        plt.show()
        print(' ')


print(' ')
"""
# TODO: may need to save the output files for each parameter set and run a cf
"""
#Measure correlation function
files = "/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/test/[01]/*/picca-gaussian-RSD-renorm-16-*.fits"
output = "test.fits.gz"
rp_max = 80. #Mpc/h
rt_max = 80. #Mpc/h
npar = 20
ntra = 20

command = '/global/homes/j/jfarr/Programs/picca/bin/do_cf.py --in-dir "dummy" --from-image {} --out {} --rp-max {} --rt-max {} --np {} --nt {} --no-project --nside 64 --nproc 64'.format(files,output,rp_max,rt_max,npar,ntra)
#process = subprocess.run(command.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)

#Run the fitter

#For now, plot the results
"""
