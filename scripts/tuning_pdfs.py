import numpy as np
import glob
import time
from multiprocessing import Pool
import multiprocessing
from astropy.io import fits
import matplotlib.pyplot as plt

from pyacolore import simulation_data, bias, utils

#base_dir = '../example_data/lya_skewers/'
base_dir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v5/v5.0.0/'
tuning_files = glob.glob('./input_files/tuning_data_a?.?_b1.65.fits')
# + glob.glob('./input_files/tuning_data_a?.?_b2.0.fits')
#tuning_files = glob.glob('./input_files/tuning_data_apow4.5_sGconst.fits')
z_values = np.array([2.0,2.2,2.4,2.6,2.8,3.0,3.2])
N_bins = 1000
z_width_value = 0.1
N_pixels = 32
f = 0.9625

plot_option = 'plot_per_z_value' #'plot_per_tuning'

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
include_thermal_effects = False
N_side = 16

def pdf_tuning(pixel_object,tuning_filename,z_values,z_width=0.2,bins=100):

    #Get tuning data
    h = fits.open(tuning_filename)
    tuning_z_values = h[1].data['z']
    tuning_alphas = h[1].data['alpha']
    tuning_betas = h[1].data['beta']
    tuning_sigma_Gs = h[1].data['sigma_G']
    n = h[1].header['n']
    k1 = h[1].header['k1']
    h.close()

    #Trim the skewers (remove low lambda cells). Exit if no QSOs are left.
    #We don't cut too tightly on the low lambda to allow for RSDs.
    lambda_buffer = 100. #A
    pixel_object.trim_skewers(lambda_min-lambda_buffer,min_catalog_z,extra_cells=1)
    if pixel_object.N_qso == 0:
        print('\nwarning: no objects left in pixel {} after trimming.'.format(pixel))
        return pixel

    #Add small scale power to the gaussian skewers:
    generator = np.random.RandomState(seed)
    def seps_z(z):
        return np.exp(np.interp(np.log(z),np.log(tuning_z_values),np.log(tuning_sigma_Gs)))
    pixel_object.add_small_scale_gaussian_fluctuations(final_cell_size,seps_z,generator,white_noise=False,lambda_min=lambda_min,IVAR_cutoff=IVAR_cutoff,n=n,k1=k1,R_kms=R_kms)

    #Remove the 'SIGMA_G' header as SIGMA_G now varies with z, so can't be stored in a header.
    sigma_G = np.sqrt(tuning_sigma_Gs**2 + measured_SIGMA_G**2)
    pixel_object.SIGMA_G = np.exp(np.interp(np.log(pixel_object.Z),np.log(tuning_z_values),np.log(sigma_G)))

    #Recompute physical skewers.
    pixel_object.compute_physical_skewers()

    #Add tau skewers to the object, starting with Lyman-alpha
    alphas = np.exp(np.interp(np.log(pixel_object.Z),np.log(tuning_z_values),np.log(tuning_alphas)))
    betas = np.exp(np.interp(np.log(pixel_object.Z),np.log(tuning_z_values),np.log(tuning_betas)))
    sigma_Gs = np.exp(np.interp(np.log(pixel_object.Z),np.log(tuning_z_values),np.log(tuning_sigma_Gs)))
    pixel_object.compute_all_tau_skewers(alphas,betas)

    #Add RSDs from the velocity skewers provided by CoLoRe.
    pixel_object.add_all_RSDs(thermal=include_thermal_effects)

    #Trim the skewers (remove low lambda cells). Exit if no QSOs are left.
    #We now cut hard at lambda min as RSDs have been implemented.
    pixel_object.trim_skewers(lambda_min,min_catalog_z,extra_cells=1)

    #Recompute alphas and betas etc to take into account trimming.
    alphas = np.exp(np.interp(np.log(pixel_object.Z),np.log(tuning_z_values),np.log(tuning_alphas)))
    betas = np.exp(np.interp(np.log(pixel_object.Z),np.log(tuning_z_values),np.log(tuning_betas)))
    sigma_Gs = np.exp(np.interp(np.log(pixel_object.Z),np.log(tuning_z_values),np.log(tuning_sigma_Gs)))

    histograms = {}
    for z_value in z_values:
        #Calculate biases.
        hist,edges = pixel_object.get_pdf_quantity('flux',z_value=z_value,z_width=z_width,bins=bins)
        z_hist = {}
        z_hist['hist'] = hist
        z_hist['edges'] = edges
        histograms[z_value] = z_hist
    
    return histograms


def pixel_tuning_bias(pixel,tuning_filename,z_values,d=0.001,z_width=0.2):

    dirname = utils.get_dir_name(base_dir,pixel)
    gaussian_filename = utils.get_file_name(dirname,'gaussian-colore',N_side,pixel)
    file_number = None
    pixel_object = simulation_data.SimulationData.get_gaussian_skewers_object(gaussian_filename,file_number,input_format,SIGMA_G=measured_SIGMA_G,IVAR_cutoff=IVAR_cutoff)
    bins = np.linspace(0.,1.,N_bins+1)

    histograms = pdf_tuning(pixel_object,tuning_filename,z_values,z_width=z_width,bins=bins)

    return histograms

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

histograms = {}
for tuning_filename in tuning_files:
    print('looking at',tuning_filename)
    tasks = [(pixel,tuning_filename,z_values,z_width_value) for pixel in pixels]

    #Run the multiprocessing pool
    if __name__ == '__main__':
        pool = Pool(processes = N_processes)
        results = []
        start_time = time.time()

        for task in tasks:
            pool.apply_async(pixel_tuning_bias,task,callback=log_result,error_callback=log_error)

        pool.close()
        pool.join()

    hist_tuning = {}
    for z_value in z_values:
        histogram = np.zeros(N_bins)
        edges = np.linspace(0.,1.,N_bins+1)
        for result in results:
            histogram += result[z_value]['hist']
        histogram /= len(results)
        hist_tuning_z_value = {}
        hist_tuning_z_value['hist'] = histogram
        hist_tuning_z_value['edges'] = edges
        hist_tuning[z_value] = hist_tuning_z_value
    histograms[tuning_filename] = hist_tuning

if plot_option == 'plot_per_z_value':
    for z_value in z_values:
        for tuning_filename in tuning_files:
            histogram = histograms[tuning_filename][z_value]['hist']
            edges = histograms[tuning_filename][z_value]['edges']
            centres = edges[:-1]/2.+edges[1:]/2.
            b_eta = np.trapz(centres*np.log(centres)*histogram,centres)
            label = tuning_filename[tuning_filename.rfind('/'):]
            plt.step(centres,histogram,where='mid',label=label+', b_eta={:2.4f}'.format(b_eta))
        plt.title(r'$z={}$'.format(z_value))
        plt.semilogy()
        plt.legend()
        plt.grid()
        plt.show()
elif plot_option == 'plot_per_tuning':
    for tuning_filename in tuning_files:
        for z_value in z_values:
            histogram = histograms[tuning_filename][z_value]['hist']
            edges = histograms[tuning_filename][z_value]['edges']
            centres = edges[:-1]/2.+edges[1:]/2.
            b_eta = np.trapz(centres*np.log(centres)*histogram,centres)
            plt.step(centres,histogram,where='mid',label=r'$z={}, b_\eta={:2.4f}$'.format(z_value,b_eta))
        title = tuning_filename[tuning_filename.rfind('/'):]
        plt.title(title)
        plt.semilogy()
        plt.legend()
        plt.grid()
        plt.show()
