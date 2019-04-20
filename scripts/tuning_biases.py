import numpy as np
import glob
import time
from multiprocessing import Pool
import multiprocessing
from astropy.io import fits

from pyacolore import simulation_data, bias, utils, tuning

#base_dir = '../example_data/lya_skewers/'
base_dir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v5/v5.0.0/'
tuning_files = glob.glob('./input_files/tuning_data_with_bias_a2.0_b1.65.fits') 
#+ glob.glob('./input_files/tuning_data_a?.?_b2.0.fits')
#tuning_files = glob.glob('./input_files/tuning_data_apow4.5_sGconst.fits')
#z_values = np.array([2.0,2.2,2.4,2.6,2.8,3.0,3.2])
z_values = np.array([2.0,2.4,2.8,3.2])
d_value = 10**-9
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
vel_boost = 1.0
include_thermal_effects = False
N_side = 16

def bias_tuning(pixel_object,tuning_filename,z_values,d=0.001,z_width=0.2,z_r0=2.5):

    #Get tuning data
    h = fits.open(tuning_filename)
    tuning_z_values = h[1].data['z']
    tuning_alphas = h[1].data['alpha']
    tuning_betas = h[1].data['beta']
    tuning_sigma_Gs = h[1].data['sigma_G']
    n = h[1].header['n']
    k1 = h[1].header['k1']
    h.close()

    transformation = tuning.transformation()
    transformation.add_parameters_from_data(tuning_z_values,tuning_alphas,tuning_betas,tuning_sigma_Gs)
    pixel_object.transformation = transformation

    #Trim the skewers (remove low lambda cells). Exit if no QSOs are left.
    #We don't cut too tightly on the low lambda to allow for RSDs.
    lambda_buffer = 100. #A
    pixel_object.trim_skewers(lambda_min-lambda_buffer,min_catalog_z,extra_cells=1)
    if pixel_object.N_qso == 0:
        print('\nwarning: no objects left in pixel {} after trimming.'.format(pixel))
        return pixel

    #Add small scale power to the gaussian skewers:
    generator = np.random.RandomState(seed)
    pixel_object.add_small_scale_gaussian_fluctuations(final_cell_size,generator,white_noise=False,lambda_min=lambda_min,IVAR_cutoff=IVAR_cutoff,n=n,k1=k1,R_kms=R_kms)
    #Recompute physical skewers.
    pixel_object.compute_physical_skewers()

    #Add tau skewers to the object, starting with Lyman-alpha
    pixel_object.compute_all_tau_skewers()

    #Add RSDs from the velocity skewers provided by CoLoRe.
    pixel_object.add_all_RSDs(thermal=include_thermal_effects)

    #Trim the skewers (remove low lambda cells). Exit if no QSOs are left.
    #We now cut hard at lambda min as RSDs have been implemented.
    pixel_object.trim_skewers(lambda_min,min_catalog_z,extra_cells=1)

    #Calculate biases.
    b = bias.get_bias_delta(pixel_object,z_values,d=d,z_width=z_width)
    b_eta = bias.get_bias_eta(pixel_object,z_values,d=d,z_width=z_width)

    return b,b_eta


def pixel_tuning_bias(pixel,tuning_filename,z_values,d=0.001,z_width=0.2,vel_boost=1.0):

    dirname = utils.get_dir_name(base_dir,pixel)
    gaussian_filename = utils.get_file_name(dirname,'gaussian-colore',N_side,pixel)
    file_number = None
    pixel_object = simulation_data.SimulationData.get_gaussian_skewers_object(gaussian_filename,file_number,input_format,SIGMA_G=measured_SIGMA_G,IVAR_cutoff=IVAR_cutoff)
    pixel_object.VEL_rows *= vel_boost

    b,b_eta = bias_tuning(pixel_object,tuning_filename,z_values,d=d,z_width=z_width)

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
    tasks = [(pixel,tuning_filename,z_values,d_value,z_width_value,vel_boost) for pixel in pixels]

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

    print(tuning_filename)
    print_string = ''
    for bs in [b_all_pixels,b_nu_all_pixels,beta_all_pixels]:
        for b in bs:
            print_string += '{:1.6f}, '.format(b)
        print_string = print_string[:-2] + '\n'
    print(print_string)
