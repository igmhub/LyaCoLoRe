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
for dir in dirs[:4]:
    pixels += [int(dir[len(dir)-dir[-2::-1].find('/')-1:-1])]


def measure_pixel_segment(pixel,z_value,alpha,beta,sigma_G_required,n,k1):
    #print('start',pixel,z_value,alpha,beta,sigma_G_required,n,k1)

    location = base_file_location + '/' + new_file_structure.format(pixel//100,pixel)

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

    ID = 0
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

def f(alpha,beta,sigma_G,n,k1):

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
        return retval

    #Define an error-tracking function.
    def log_error(retval):
        print('Error:',retval)

    ################################################################################


    print('looking at params: a={:2.4f}, b={:2.4f}, sG={:2.4f}, n={:2.4f}, k1={:2.6f}'.format(alpha,beta,sigma_G,n,k1))

    tasks = [(pixel,z_value,alpha,beta,sigma_G,n,k1) for pixel in pixels for z_value in z_values]

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

    overall_chi2 = 0.
    for m in combined_pixels_set.measurements:
        m.add_mean_F_chi2(eps=0.05)
        m.add_Pk1D_chi2(max_k=max_k)
        m.add_total_chi2()
        overall_chi2 += m.total_chi2

    """
    combined_z_pixels_set = measurement_set.combine_zs()

    if len(combined_z_pixels_set.measurements) == 1:
        m = combined_z_pixels_set.measurements[0]
    """
    print('overall chi2={:2.4f}'.format(overall_chi2))
    print(' ')
    return overall_chi2

m = Minuit(f,alpha=0.82,error_alpha=0.05,limit_alpha=(0.,10.),
             beta=1.65,error_beta=0.05,limit_beta=(0.,10.),fix_beta=True,
             sigma_G=4.83,error_sigma_G=0.05,limit_sigma_G=(0.,20.),
             n=0.7,error_n=0.05,limit_n=(0.,10.),fix_n=False,
             k1=0.001,error_k1=0.00005,limit_k1=(0.,0.1),fix_k1=False,
             )

m.print_param()
m.migrad()
m.print_param()
