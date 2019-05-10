import numpy as np
import glob
import time
from multiprocessing import Pool
import multiprocessing
from astropy.io import fits

from lyacolore import simulation_data, bias, tuning, utils

#base_dir = '../example_data/lya_skewers/'
base_dir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v5/v5.0.0/'
tuning_files = glob.glob('./input_files/tuning_data_a?.?_b1.65.fits')
#+ glob.glob('./input_files/tuning_data_a?.?_b2.0.fits')
#tuning_files = glob.glob('./input_files/tuning_data_apow4.5_sGconst.fits')
z_values = np.array([2.0,2.2,2.4,2.6,2.8,3.0,3.2])
d_value = 10**-3
z_width_value = 0.1
N_pixels = 1
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
include_thermal_effects = False
N_side = 16

pixel = 0
dirname = utils.get_dir_name(base_dir,pixel)
gaussian_filename = utils.get_file_name(dirname,'gaussian-colore',N_side,pixel)
file_number = None
pixel_object = simulation_data.SimulationData.get_gaussian_skewers_object(gaussian_filename,file_number,input_format,SIGMA_G=measured_SIGMA_G,IVAR_cutoff=IVAR_cutoff)

tuning_filename = tuning_files[0]
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

#Add small scale power to the gaussian skewers:
generator = np.random.RandomState(seed)
pixel_object.add_small_scale_gaussian_fluctuations(final_cell_size,generator,white_noise=False,lambda_min=lambda_min,IVAR_cutoff=IVAR_cutoff,n=n,k1=k1,R_kms=R_kms)

#Recompute physical skewers.
pixel_object.compute_physical_skewers()

#Add tau skewers to the object, starting with Lyman-alpha
pixel_object.compute_all_tau_skewers()

import numpy as np
import copy

from pyacolore import utils

lya = utils.lya_rest
d = 10**-3

data = copy.deepcopy(pixel_object)
data_noRSDs = copy.deepcopy(data)

data.add_all_RSDs()

z_value = 2.4
z_width = 0.2

z_min = z_value - 0.5*z_width
z_max = z_value + 0.5*z_width
lambda_min = lya * (1 + z_min)
lambda_max = lya * (1 + z_max)

lambda_buffer = 100. #A
min_catalog_z = 1.8

data_noRSDs_z_val = copy.deepcopy(data_noRSDs)
data_noRSDs_z_val.trim_skewers(lambda_min-lambda_buffer,min_catalog_z,lambda_max=lambda_max+lambda_buffer,extra_cells=1)

grad_increase = copy.deepcopy(data_noRSDs_z_val)
grad_increase.add_all_RSDs(thermal=include_thermal_effects,d=d,z_r0=z_value)
grad_increase.trim_skewers(lambda_min,min_catalog_z,lambda_max=lambda_max)

grad_decrease = copy.deepcopy(data_noRSDs_z_val)
grad_decrease.add_all_RSDs(thermal=include_thermal_effects,d=-d,z_r0=z_value)
grad_decrease.trim_skewers(lambda_min,min_catalog_z,lambda_max=lambda_max)

data_z_val = copy.deepcopy(data)
data_z_val.trim_skewers(lambda_min,min_catalog_z,lambda_max=lambda_max)

data_noRSDs_z_val.trim_skewers(lambda_min,min_catalog_z,lambda_max=lambda_max)

#We get means across the z-chunk and combine once bias has been computed.
#This avoids overweighting the low end of the chunk.
mean_F_increase = grad_increase.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)
mean_F_decrease = grad_decrease.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)
mean_F = data_z_val.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)

bias = np.average((1/mean_F) * (1/(2.*d)) * (mean_F_increase - mean_F_decrease))

####################################

plt.plot(grad_increase.Z,grad_increase.lya_absorber.transmission()[2,:],label='increase')
plt.plot(data_z_val.Z,data_z_val.lya_absorber.transmission()[2,:],label='original')
plt.plot(grad_decrease.Z,grad_decrease.lya_absorber.transmission()[2,:],label='decrease')
plt.xlim(2.30,2.32)
plt.legend();plt.grid();plt.show()

for key in grad_increase.lya_absorber.__dict__.keys():
    print(key)
    grad_increase.lya_absorber.__dict__[key] == grad_decrease.lya_absorber.__dict__[key]
    print(' ')

plt.plot(grad_increase.Z,mean_F_increase,label='increase')
plt.plot(data_z_val.Z,mean_F,label='original')
plt.plot(grad_decrease.Z,mean_F_decrease,label='decrease')
plt.xlim(2.30,2.32)
plt.legend();plt.grid();plt.show()

