import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
import subprocess
from multiprocessing import Pool
import multiprocessing
import time

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
z_values = [2.5,3.0,3.5]
z_width = 0.2

cell_size = 0.25 #Mpc/h

#Open up the Gaussian colore files
base_file_location = '/Users/jfarr/Projects/test_data/test/'
N_side = 16

new_file_structure = '{}/{}/'               #pixel number//100, pixel number
new_filename_structure = '{}-{}-{}.fits'    #file type, nside, pixel number

input_format = 'gaussian_colore'

# TODO: get pixels from those directories created by make_master.py
pixels = list(range(100))
pixels=[0]

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

beta = 1.65
betas = beta*np.ones(tuning_z_values.shape[0])

#Get the values of alpha, beta, sigma_G for the z_values
alpha_values = []
beta_values = []
sigma_G_values = []

for z_value in z_values:
    alpha_values += [np.interp(z_value,tuning_z_values,alphas)]
    beta_values += [np.interp(z_value,tuning_z_values,betas)]
    sigma_G_values += [np.interp(z_value,tuning_z_values,desired_sigma_G_values)]

multipliers = np.linspace(0.5,1.5,5)

#Extract the values of parameters to optimies over
parameters_list = []

lookup = {}

import itertools
for i,z_value in enumerate(z_values):
    a = alpha_values[i] * multipliers
    b = beta_values[i] * multipliers
    sG = sigma_G_values[i] * multipliers
    parameters_list += list(itertools.product([z_value],a,b,sG))

    parameter_values_list = list(itertools.product(a,b,sG))

    # TODO: need a better ID number of some kind
    ID = 0
    z_dict = {}
    for parameter_values in parameter_values_list:
        parameters_dict = {'alpha':parameter_values[0],'beta':parameter_values[1],'sigma_G':parameter_values[2]}
        results_dict = {'mean_F':None,'A_F':None,'B_F':None}
        ID_dict = {'parameters':parameters_dict,'results':results_dict}
        z_dict = {**z_dict,**{ID:ID_dict}}
        ID += 1


    data_type = [('alpha','f4'),('beta','f4'),('sigma_G','f4')]
    parameters = np.array(parameter_values,dtype=data_type)

    lookup = {**lookup,**{z_value:z_dict}}


print(lookup)

def measure_pixel_segment(pixel,z_value,alpha,beta,sigma_G_required):

    input = (pixel,z_value,alpha,beta,sigma_G_required)

    location = base_file_location + '/' + new_file_structure.format(pixel//100,pixel)
    mean_F_data = np.array(list(zip(tuning_z_values,desired_mean_F)))

    # TODO: this is a problem
    measured_SIGMA_G = 1.17

    #We work from the gaussian colore files made in 'pixelise gaussian skewers'.
    gaussian_filename = new_filename_structure.format('gaussian-colore',N_side,pixel)

    #Make a pixel object from it.
    data = pixelise.simulation_data.get_gaussian_skewers_object(location+gaussian_filename,None,input_format,SIGMA_G=measured_SIGMA_G,IVAR_cutoff=IVAR_cutoff)

    #Determine the desired sigma_G by sampling
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
    data.add_small_scale_gaussian_fluctuations(cell_size,data.Z,np.ones(data.Z.shape[0])*extra_sigma_G,generator,amplitude=1.0,white_noise=False,lambda_min=0.0,IVAR_cutoff=lya)

    #Convert to flux
    data.compute_physical_skewers()
    data.compute_tau_skewers(alpha=np.ones(data.Z.shape[0])*alpha,beta=beta)
    data.add_RSDs(np.ones(data.Z.shape[0])*alpha,beta,thermal=False)
    data.compute_flux_skewers()

    #Trim the skewers again to get rid of the additional cells
    lambda_min_val = lya*(1 + z_value - z_width/2)
    lambda_max_val = lya*(1 + z_value + z_width/2)
    data.trim_skewers(lambda_min_val,min_cat_z,lambda_max=lambda_max_val,whole_lambda_range=True)

    #Measure mean flux
    mean_F = data.get_mean_flux(z_value=z_value)

    #Fit model
    mean_F_model = tuning.get_mean_F_model(z_value)

    #Measure P1D
    #Do i want this to be a pixelise function?
    #this is returning complex results atm
    mean_F = np.average(data.F_rows)
    delta_F_rows = data.F_rows/mean_F - 1

    k_kms, Pk_kms, var_kms = Pk1D.get_Pk1D(delta_F_rows,data.IVAR_rows,data.R,data.Z,z_value,z_width=0.2,N_processes=1)

    Pk1D_results = (k_kms,Pk_kms,var_kms)

    #Fit model
    def model_Pk_kms(k_kms,A_F,B_F):
        return tuning.P1D_z_kms_PD2013(k_kms,z_value,A_F=A_F,B_F=B_F)
    fit = curve_fit(model_Pk_kms,k_kms,Pk_kms,p0=(0.064,3.55))
    #something to do with comparing to the default values

    model_Pk_kms_fit = model_Pk_kms(k_kms,fit[0][0],fit[0][1])
    model_Pk_kms_default = model_Pk_kms(k_kms,0.064,3.55)

    Pk1D_results_model = (k_kms,fit,model_Pk_kms_fit,model_Pk_kms_default)

    #Convert back to small cells
    #need a new function to merge cells back together

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


    cf_location = ''
    result = (mean_F,mean_F_model,Pk1D_results,Pk1D_results_model,cf_location)

    return (input,result)

tasks = [(pixel,*parameters) for parameters in parameters_list for pixel in pixels]

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = N_processes)
    start_time = time.time()
    results = []

    for task in tasks:
        pool.apply_async(measure_pixel_segment,task,callback=log_result,error_callback=log_error)

    pool.close()
    pool.join()

# TODO: write something to combine the pixels together
results_dict = {}
for result in results:
    results_dict = {**results_dict,**{result[0]:result[1]}}

for z_value in z_values:
    # TODO: need a way of labelling rather than using a numerical index
    # TODO: need some kind of measurement ID to indicate that the same parameters at the same z value are being used. Maybe need a lookup at the start, so then don't need to pass parameter values around so much. Then I can just have pixel,parameters and take all pixels that match a certain parameter
    relevant_keys = [key for key in results_dict.keys() if key[1] == z_value]





errors = []

for i,result in enumerate(results):
    z_value = results[i][0][1]
    alpha = results[i][0][2]
    beta = results[i][0][3]
    sigma_G_required = results[i][0][4]

    mean_F = results[i][1][0]
    mean_F_model = results[i][1][1]

    k_kms = results[i][1][2][0]
    Pk_kms = results[i][1][2][1]
    var_kms = results[i][1][2][2]

    k_kms = results[i][1][3][0]
    fit = results[i][1][3][1]
    model_Pk_kms_fit = results[i][1][3][2]
    model_Pk_kms_default = results[i][1][3][3]

    mean_F_error = mean_F/mean_F_model - 1
    A_F_error = fit[0][0]/0.064 - 1
    B_F_error = fit[0][1]/3.55 - 1

    error = np.sqrt(mean_F_error**2 + A_F_error**2 + B_F_error**2)
    errors += [error]

best_result = np.argmin(errors)

z_value = results[best_result][0][1]
alpha = results[best_result][0][2]
beta = results[best_result][0][3]
sigma_G_required = results[best_result][0][4]

mean_F = results[best_result][1][0]
mean_F_model = results[best_result][1][1]

k_kms = results[best_result][1][2][0]
Pk_kms = results[best_result][1][2][1]
var_kms = results[best_result][1][2][2]

k_kms = results[best_result][1][3][0]
fit = results[best_result][1][3][1]
model_Pk_kms_fit = results[best_result][1][3][2]
model_Pk_kms_default = results[best_result][1][3][3]

mean_F_error = mean_F/mean_F_model - 1
A_F_error = fit[0][0]/0.064 - 1
B_F_error = fit[0][1]/3.55 - 1
print('mean F error is {:3.0%}'.format(mean_F_error))
print('A_F error is {:3.0%}'.format(A_F_error))
print('B_F error is {:3.0%}'.format(B_F_error))

plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.title('alpha={:2.2f}, beta={:2.2f}, sG={:2.2f}'.format(alpha,beta,sigma_G_required))
plt.errorbar(k_kms,Pk_kms,yerr=np.sqrt(var_kms),fmt='o',label='measured',color='orange')
plt.plot(k_kms,model_Pk_kms_fit,label='model fit: A_F={:2.2f}, B_F={:2.2f}'.format(fit[0][0],fit[0][1]),color='b')
plt.plot(k_kms,model_Pk_kms_default,label='model default: A_F={:2.2f}, B_F={:2.2f}'.format(0.064,3.55),color=(0.5,0.5,0.5))
plt.fill_between(k_kms,0.9*model_Pk_kms_default,1.1*model_Pk_kms_default,label='model default +/- 10%',color=(0.5,0.5,0.5),alpha=0.5)
#    plt.plot(k_kms,independent.power_kms(z_value,k_kms,cell_size*general.get_dkms_dhMpc(z_value),False),label='added')
plt.semilogy()
plt.semilogx()
plt.ylabel('Pk1D')
plt.xlabel('k / kms-1')
plt.legend()
plt.grid()
plt.savefig('Pk1D_abs_slope1.5_alpha{:2.2f}_beta{:2.2f}_sG{:2.2f}.pdf'.format(alpha,beta,sigma_G_required))
plt.show()



"""
mean_Pk_kms = np.average([result[1] for result in results],axis=0)

def model_Pk_kms(k_kms,A_F,B_F):
    return tuning.P1D_z_kms_PD2013(k_kms,2.5,A_F=A_F,B_F=B_F)
fit = curve_fit(model_Pk_kms,k_kms,mean_Pk_kms,p0=(0.064,3.55))
#something to do with comparing to the default values

model_Pk_kms_fit = model_Pk_kms(k_kms,fit[0][0],fit[0][1])
model_Pk_kms_default = model_Pk_kms(k_kms,0.064,3.55)


plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.errorbar(k_kms,mean_Pk_kms,fmt='o',label='measured',color='orange')
plt.plot(k_kms,model_Pk_kms_fit,label='model fit: A_F={:2.2f}, B_F={:2.2f}'.format(fit[0][0],fit[0][1]),color='b')
plt.plot(k_kms,model_Pk_kms_default,label='model default: A_F={:2.2f}, B_F={:2.2f}'.format(0.064,3.55),color=(0.5,0.5,0.5))
plt.fill_between(k_kms,0.9*model_Pk_kms_default,1.1*model_Pk_kms_default,label='model default +/- 10%',color=(0.5,0.5,0.5),alpha=0.5)
plt.plot(k_kms,independent.power_kms(2.5,k_kms,cell_size*general.get_dkms_dhMpc(2.5),False),label='added')
plt.semilogy()
plt.semilogx()
plt.legend()
plt.grid()
plt.savefig('Pk1D_abs_slope1.5.pdf')
plt.show()
"""
