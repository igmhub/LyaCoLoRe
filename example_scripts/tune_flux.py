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
for dir in dirs[:3]:
    pixels += [int(dir[len(dir)-dir[-2::-1].find('/')-1:-1])]
#pixels=[0]

################################################################################

"""
Define the multiprocessing tracking functions
"""

#Define a progress-tracking function.
def log_result(retval):

    results.append(retval[0])
    new_results.append(retval[1])
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

multipliers = np.linspace(0.8,1.2,2)

#Extract the values of parameters to optimise over
parameters_list = []

lookup = {}

import itertools
for i,z_value in enumerate(z_values):
    a = [alpha_values[i]] * multipliers
    b = [beta_values[i]]# * multipliers
    sG = [sigma_G_values[i]] * multipliers

    n = [n_values[i]] * multipliers
    k1 = [k1_values[i]] * multipliers

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

    #Measure P1D
    mean_F = np.average(data.F_rows)
    delta_F_rows = data.F_rows/mean_F - 1

    #Do i want this to be a pixelise function?
    k_kms, Pk_kms, var_kms = Pk1D.get_Pk1D(delta_F_rows,data.IVAR_rows,data.R,data.Z,z_value,z_width=0.2,N_processes=1)

    Pk1D_results = (k_kms,Pk_kms,var_kms)

    #Measure mean flux
    mean_F = data.get_mean_flux(z_value=z_value,z_width=z_width)

    #Convert back to small cells for cf measurement
    #need a new function to merge cells back together

    # TODO: this doesn't want to go here - we need to save the files and then measure it.
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

    data_type = [('z_value','f4'),('ID','i4'),('pixel','i4'),('N','i4'),('mean_F','f4'),('k_kms','O'),('Pk_kms','O')]
    result = np.array([(z_value,ID,pixel,data.N_qso,mean_F,k_kms,Pk_kms)],dtype=data_type)

    measurement = tuning.measurement(ID,z_value,z_width,data.N_qso,n,k1,alpha,beta,sigma_G_required,pixels=[pixel])
    measurement.add_mean_F_measurement(data)
    measurement.add_Pk1D_measurement(data)
    #measurement_set.add_measurement(measurement)

    return (result,measurement)


tasks = [(pixel,z_value,ID,lookup) for ID in ID_list for pixel in pixels for z_value in z_values]

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = N_processes)
    start_time = time.time()
    results = []
    new_results = []

    for task in tasks:
        pool.apply_async(measure_pixel_segment,task,callback=log_result,error_callback=log_error)

    pool.close()
    pool.join()


measurement_set = tuning.measurement_set(measurements=new_results)


print('originally number measuremts is:',len(measurement_set.measurements))
combined_set = measurement_set.combine_pixels()
print('after combining, number measuremts is:',len(combined_set.measurements))

#Fit model
def get_model_Pk_kms(k_kms,A_F,B_F):
    return tuning.P1D_z_kms_PD2013(k_kms,z_value,A_F=A_F,B_F=B_F)


for z_value in z_values:
    print('z value: {}'.format(z_value))
    z_set = measurement_set.z_filter(z_value)
    for s_parameter_values in s_parameter_values_list:
        n = s_parameter_values[0]
        k1 = s_parameter_values[1]
        z_s_set = z_set.s_filter(n,k1)
        for m in z_s_set.measurements:
            m.add_Pk1D_error(max_k=max_k)
            m.add_mean_F_error()
            m.add_Pk1D_chi2(max_k=max_k)
            m.add_mean_F_chi2()
        best = z_s_set.get_best_measurement()
        max_j = np.searchsorted(best.k_kms,max_k)
        fit = curve_fit(get_model_Pk_kms,best.k_kms[:max_j],best.Pk_kms[:max_j],p0=(0.064,3.55))

        print('alpha = {:2.2f}'.format(best.alpha))
        print('beta = {:2.2f}'.format(best.beta))
        print('sigma_G = {:2.2f}'.format(best.sigma_G))
        print('n = {:2.2f}'.format(best.n))
        print('k1 = {:2.4f}'.format(best.k1))
        print(' ')
        print('mean F = {:2.2f}, error is {:3.0%}'.format(best.mean_F,best.mean_F_error))
        print('Pk1D average error is {:3.0%}'.format(best.Pk_kms_error))
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


print('############################')

results = np.array(results)

mean_F_errors = np.array((multipliers.shape[0],multipliers.shape[0]))
Pk1D_errors = np.array((multipliers.shape[0],multipliers.shape[0]))
A_F_errors = np.array((multipliers.shape[0],multipliers.shape[0]))
B_F_errors = np.array((multipliers.shape[0],multipliers.shape[0]))

#s_errors_grid = np.zeros()
#t_errors_grid

for z_value in z_values:
    min_error = 1.0
    z_results = results[results['z_value']==z_value]
    print('z value: {}'.format(z_value))

    for ID in ID_list:

        z_ID_results = z_results[z_results['ID']==ID]

        lookup[z_value][ID]['results']['mean_F'] = np.average(z_ID_results['mean_F'],weights=z_ID_results['N'])

        # TODO: need to check the k_kms are all the same
        lookup[z_value][ID]['results']['Pk1D']['k_kms'] = z_ID_results['k_kms'][0]

        # TODO: these ought really to be weighted. Doesn't make much differences as the pixels all have the same number of skewers roughly, but it may be relevant later on
        lookup[z_value][ID]['results']['Pk1D']['Pk_kms'] = np.average(z_ID_results['Pk_kms'],axis=0)
        lookup[z_value][ID]['results']['Pk1D']['var_kms'] = np.average((z_ID_results['Pk_kms']**2),axis=0) - lookup[z_value][ID]['results']['Pk1D']['Pk_kms']**2

        #Need to measure cf here
        lookup[z_value][ID]['results']['bias'] = None
        lookup[z_value][ID]['results']['beta'] = None

        #Get model values
        max_j = np.searchsorted(lookup[z_value][ID]['results']['Pk1D']['k_kms'],max_k)
        fit = curve_fit(get_model_Pk_kms,lookup[z_value][ID]['results']['Pk1D']['k_kms'][:max_j],lookup[z_value][ID]['results']['Pk1D']['Pk_kms'][:max_j],p0=(0.064,3.55))
        A_F = fit[0][0]
        B_F = fit[0][1]
        lookup[z_value][ID]['results']['A_F'] = A_F
        lookup[z_value][ID]['results']['B_F'] = B_F

        # TODO: some sort of goodness of fit estimator may be better here?
        A_F = 0.064
        B_F = 3.55
        #model_Pk_kms = get_model_Pk_kms(lookup[z_value][ID]['results']['Pk1D']['k_kms'],A_F,B_F)
        model_Pk_kms = get_model_Pk_kms(lookup[z_value][ID]['results']['Pk1D']['k_kms'],0.064,3.55)
        Pk_differences = lookup[z_value][ID]['results']['Pk1D']['Pk_kms'][:max_j] - model_Pk_kms[:max_j]
        Pk_error = np.average(abs(Pk_differences)/model_Pk_kms[:max_j])

        mean_F_model = tuning.get_mean_F_model(z_value)

        lookup[z_value][ID]['errors']['mean_F'] = lookup[z_value][ID]['results']['mean_F']/mean_F_model - 1
        lookup[z_value][ID]['errors']['Pk1D'] = Pk_error
        lookup[z_value][ID]['errors']['A_F'] = lookup[z_value][ID]['results']['A_F']/0.064 - 1
        lookup[z_value][ID]['errors']['B_F'] = lookup[z_value][ID]['results']['B_F']/3.55 - 1
        lookup[z_value][ID]['errors']['bias'] = 0.
        lookup[z_value][ID]['errors']['beta'] = 0.

        error = np.sqrt(np.sum([lookup[z_value][ID]['errors'][key]**2 for key in lookup[z_value][ID]['errors'].keys()]))

        #print(lookup[z_value][ID]['parameters'])
        #print(error)
        #print(' ')

        if error < min_error:
            min_error = error
            min_err_ID = ID

    best = lookup[z_value][min_err_ID]

    print('Best fit is:')
    print('alpha = {:2.2f}'.format(best['t_parameters']['alpha']))
    print('beta = {:2.2f}'.format(best['t_parameters']['beta']))
    print('sigma_G = {:2.2f}'.format(best['t_parameters']['sigma_G']))
    print('n = {:2.2f}'.format(best['s_parameters']['n']))
    print('k1 = {:2.4f}'.format(best['s_parameters']['k1']))
    print(' ')
    print('mean F = {:2.2f}, error is {:3.0%}'.format(best['results']['mean_F'],best['errors']['mean_F']))
    print('Pk1D average error is {:3.0%}'.format(best['errors']['Pk1D']))
    print('A_F = {:2.2f}, error is {:3.0%}'.format(best['results']['A_F'],best['errors']['A_F']))
    print('B_F = {:2.2f}, error is {:3.0%}'.format(best['results']['B_F'],best['errors']['B_F']))

    plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
    plt.title('z_value={:2.2f}: alpha={:2.2f}, beta={:2.2f}, sG={:2.2f}, n={:2.2f}, k1={:2.4f}'.format(z_value,best['t_parameters']['alpha'],best['t_parameters']['beta'],best['t_parameters']['sigma_G'],best['s_parameters']['n'],best['s_parameters']['k1']))
    plt.errorbar(best['results']['Pk1D']['k_kms'],best['results']['Pk1D']['Pk_kms'],yerr=np.sqrt(best['results']['Pk1D']['var_kms']),fmt='o',label='measured',color='orange')
    plt.plot(best['results']['Pk1D']['k_kms'],get_model_Pk_kms(best['results']['Pk1D']['k_kms'],best['results']['A_F'],best['results']['B_F']),label='model fit: A_F={:2.2f}, B_F={:2.2f}'.format(best['results']['A_F'],best['results']['B_F']),color='b')
    plt.plot(best['results']['Pk1D']['k_kms'],get_model_Pk_kms(best['results']['Pk1D']['k_kms'],0.064,3.55),label='model default: A_F={:2.2f}, B_F={:2.2f}'.format(0.064,3.55),color=(0.5,0.5,0.5))
    plt.axvline(x=max_k,c=(0.,0.,0.),label='max fitting k value')
    #plt.fill_between(k_kms,0.9*model_Pk_kms_default,1.1*model_Pk_kms_default,label='model default +/- 10%',color=(0.5,0.5,0.5),alpha=0.5)
    #plt.plot(k_kms,independent.power_kms(z_value,k_kms,cell_size*general.get_dkms_dhMpc(z_value),False),label='added')
    plt.semilogy()
    plt.semilogx()
    plt.ylabel('Pk1D')
    plt.xlabel('k / kms-1')
    plt.legend()
    plt.grid()
    plt.savefig('test.pdf')
    #plt.savefig('Pk1D_abs_slope1.5_alpha{:2.2f}_beta{:2.2f}_sG{:2.2f}.pdf'.format(alpha,beta,sigma_G_required))
    plt.show()
    print(' ')




"""
for z_value in z_values:
    # TODO: need a way of labelling rather than using a numerical index
    # TODO: need some kind of measurement ID to indicate that the same parameters at the same z value are being used. Maybe need a lookup at the start, so then don't need to pass parameter values around so much. Then I can just have pixel,parameters and take all pixels that match a certain parameter
    relevant_keys = [key for key in results_dict.keys() if key[1] == z_value]
"""

"""
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
"""

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
