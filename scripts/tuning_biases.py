import numpy as np
import glob
import time
from multiprocessing import Pool
import multiprocessing
from astropy.io import fits
import matplotlib.pyplot as plt

from lyacolore import simulation_data, bias, utils, tuning

lya = utils.lya_rest

#base_dir = '../example_data/lya_skewers/'
base_dir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v9/v9.0.0_all_files/'
tuning_files = glob.glob('./input_files/tuning_data_with_bias_vel1.3_b1.65_lr1200.fits')
#tuning_files = glob.glob('./input_files/tuning_data_with_bias_newRSD_vel1.0_afree_b1.65.fits')
#+ glob.glob('./input_files/tuning_data_a?.?_b2.0.fits')
#tuning_files = glob.glob('./input_files/tuning_data_apow4.5_sGconst.fits')
#z_values = np.array([2.0,2.2,2.4,2.6,2.8,3.0,3.2])
#z_values = np.array([2.0,2.4,2.8,3.2])
z_values = np.array([2.4])
d_values = 10**np.linspace(-11,-1,6)
z_width_value = 0.1
N_pixels = 32
f = 0.9625
z_r0 = 2.5

#pixels = np.random.choice(3072,size=N_pixels,replace=False)
pixels = list(range(N_pixels))
N_processes = np.min((N_pixels,64))

#Define parameters.
seed = 123
input_format = 'gaussian_colore'
measured_SIGMA_G = 1.165
IVAR_cutoff = 1200. #A
lambda_min = 3550.
min_catalog_z = 1.8
final_cell_size = 0.25
R_kms = 25.
a_v = 1.3
include_thermal_effects = False
N_side = 16

def bias_tuning(pixel_object,tuning_filename,z_values,d_delta=10**-9,d_eta=10**-2,z_width=0.2,z_r0=2.5):

    #Get tuning data
    h = fits.open(tuning_filename)
    tuning_z_values = h[1].data['z']
    tuning_alphas = h[1].data['alpha']
    tuning_betas = h[1].data['beta']
    tuning_sigma_Gs = h[1].data['sigma_G']
    n = h[1].header['n']
    k1 = h[1].header['k1']
    h.close()

    transformation = tuning.Transformation()
    transformation.add_parameters_from_data(tuning_z_values,tuning_alphas,tuning_betas,tuning_sigma_Gs,n,k1,R_kms,a_v)
    pixel_object.transformation = transformation

    pixel_object.scale_velocities(use_transformation=True)

    #Trim the skewers (remove low lambda cells). Exit if no QSOs are left.
    #We don't cut too tightly on the low lambda to allow for RSDs.
    lambda_buffer = 100. #A
    pixel_object.trim_skewers(lambda_min-lambda_buffer,min_catalog_z,extra_cells=1)
    if pixel_object.N_qso == 0:
        print('\nwarning: no objects left in pixel {} after trimming.'.format(pixel))
        return pixel

    #trim skewers to the minimal length
    #extra = 4.0
    #z_lower_cut = np.min(z_values) - z_width*(1+extra)/2.
    #z_upper_cut = np.max(z_values) + z_width*(1+extra)/2.
    #lambda_min_val = np.min([lambda_min,lya*(1 + z_lower_cut)])
    #lambda_max_val = lya*(1 + z_upper_cut)
    #pixel_object.trim_skewers(lambda_min_val,min_catalog_z,lambda_max=lambda_max_val,whole_lambda_range=False)

    #Add small scale power to the gaussian skewers:
    generator = np.random.RandomState(seed)
    pixel_object.add_small_scale_fluctuations(final_cell_size,generator,white_noise=False,lambda_min=0.0,IVAR_cutoff=IVAR_cutoff,use_transformation=True)
    #Recompute physical skewers.
    pixel_object.compute_physical_skewers()

    #Add tau skewers to the object, starting with Lyman-alpha
    pixel_object.compute_all_tau_skewers()

    #Get RSD weights.
    pixel_object.compute_RSD_weights()

    #Add RSDs from the velocity skewers provided by CoLoRe.
    pixel_object.add_all_RSDs(thermal=include_thermal_effects)

    #Trim the skewers (remove low lambda cells). Exit if no QSOs are left.
    #We now cut hard at lambda min as RSDs have been implemented.
    #pixel_object.trim_skewers(lambda_min,min_catalog_z,extra_cells=1)

    #Calculate biases.
    b = bias.get_bias_delta(pixel_object,z_values,weights=pixel_object.RSD_weights,d=d_delta,z_width=z_width)
    b_eta = bias.get_bias_eta(pixel_object,z_values,d=d_eta,z_width=z_width)

    return b,b_eta


def pixel_tuning_bias(pixel,tuning_filename,z_values,d_delta=10**-9,d_eta=10**-2,z_width=0.2,a_v=1.0):

    dirname = utils.get_dir_name(base_dir,pixel)
    gaussian_filename = utils.get_file_name(dirname,'gaussian-colore',N_side,pixel,compressed=True)
    file_number = None
    pixel_object = simulation_data.SimulationData.get_skewers_object(gaussian_filename,file_number,input_format,SIGMA_G=measured_SIGMA_G,IVAR_cutoff=IVAR_cutoff)

    b,b_eta = bias_tuning(pixel_object,tuning_filename,z_values,d_delta=d_delta,d_eta=d_eta,z_width=z_width)

    return (b,b_eta)

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

for tuning_filename in tuning_files:
    print(tuning_filename)

    final_bs = np.zeros((z_values.shape[0],3,d_values.shape[0]))

    for i,d_value in enumerate(d_values):
        print(' ->',d_value)

        d_delta = d_value
        d_eta = d_value

        tasks = [(pixel,tuning_filename,z_values,d_delta,d_eta,z_width_value,a_v) for pixel in pixels]

        #Run the multiprocessing pool
        if __name__ == '__main__':
            pool = Pool(processes = N_processes)
            results = []
            start_time = time.time()

            for task in tasks:
                pool.apply_async(pixel_tuning_bias,task,callback=log_result,error_callback=log_error)

            pool.close()
            pool.join()

        biases_all_pixels = np.array(results)

        b_all_pixels = np.average(biases_all_pixels[:,0,:],axis=0)
        b_nu_all_pixels = np.average(biases_all_pixels[:,1,:],axis=0)
        beta_all_pixels = b_nu_all_pixels * f / b_all_pixels

        final_bs[:,0,i] = b_all_pixels
        final_bs[:,1,i] = b_nu_all_pixels
        final_bs[:,2,i] = beta_all_pixels

        print_string = ''
        for bs in [b_all_pixels,b_nu_all_pixels,beta_all_pixels]:
            for b in bs:
                print_string += '{:1.6f}, '.format(b)
            print_string = print_string[:-2] + '\n'
        print(print_string)

    for i,z_value in enumerate(z_values):
        plt.plot(d_values,final_bs[i,0,:],label=r'$b_\delta$')
        plt.plot(d_values,final_bs[i,1,:],label=r'$b_\eta$')
        #plt.plot(d_values,final_bs[i,2,:],label=r'$b_\delta$')
        plt.title(tuning_filename+', z ='+str(z_value))
        plt.legend()
        plt.grid()
        plt.semilogx()
        plt.show()
