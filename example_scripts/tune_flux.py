import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
import subprocess

import pixelise
import Pk1D
import tuning
import independent
import general

lya = 1215.67

lambda_min = 3550.0
min_cat_z = 1.8
IVAR_cutoff = 1150.0

#Get the starting values of alpha, beta and sigma_G from file
#Decide which z values we are going to tune
z_values = [2.5]
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

tasks = [(pixel) for pixel in pixels]

def measure_pixel(pixel):

    #load initial tuning data
    tuning_file = 'input_files/tune_small_scale_fluctuations.fits'
    h = fits.open(tuning_file)
    tune_small_scale_fluctuations = h['DATA'].data
    h.close()
    tuning_z_values = tune_small_scale_fluctuations['z']
    desired_sigma_G_values = tune_small_scale_fluctuations['sigma_G']
    desired_mean_F = tune_small_scale_fluctuations['mean_F']
    alphas = tune_small_scale_fluctuations['alpha']
    beta = 2.0*1.65

    location = base_file_location + '/' + new_file_structure.format(pixel//100,pixel)
    mean_F_data = np.array(list(zip(tuning_z_values,desired_mean_F)))

    # TODO: this is a problem
    measured_SIGMA_G = 1.17

    #We work from the gaussian colore files made in 'pixelise gaussian skewers'.
    gaussian_filename = new_filename_structure.format('gaussian-colore',N_side,pixel)

    #Make a pixel object from it.
    data = pixelise.simulation_data.get_gaussian_skewers_object(location+gaussian_filename,None,input_format,SIGMA_G=measured_SIGMA_G,IVAR_cutoff=IVAR_cutoff)

    #Determine the desired sigma_G by sampling
    extra_sigma_G_values = np.sqrt(desired_sigma_G_values**2 - measured_SIGMA_G**2)

    #trim skewers
    data.trim_skewers(lambda_min,min_cat_z,extra_cells=1)

    #add small scale fluctuations
    seed = int(str(N_side) + str(pixel))
    generator = np.random.RandomState(seed)
    data.add_small_scale_gaussian_fluctuations(cell_size,tuning_z_values,extra_sigma_G_values,generator,amplitude=1.0,white_noise=False,lambda_min=0.0,IVAR_cutoff=lya)

    #Convert to flux
    data.compute_physical_skewers()
    data.compute_tau_skewers(alpha=np.interp(data.Z,tuning_z_values,alphas),beta=beta)
    data.add_RSDs(np.interp(data.Z,tuning_z_values,alphas),beta,thermal=False)
    data.compute_flux_skewers()

    mean_F_results = []
    Pk1D_results = []
    cf_results = []

    #For each z value, measure mean Flux, Pk1D and correlation function
    for z_value in z_values:
        #duplicate the data and crop it
        cropped_data = data
        lambda_min_val = lya*(1 + z_value - z_width/2)
        lambda_max_val = lya*(1 + z_value + z_width/2)
        cropped_data.trim_skewers(lambda_min_val,min_cat_z,lambda_max=lambda_max_val,whole_lambda_range=True)

        #Measure mean flux
        #make function to measure mean flux at a given z

        #Fit model

        #Measure P1D
        #Do i want this to be a pixelise function?
        #this is returning complex results atm
        mean_F = np.average(cropped_data.F_rows)
        delta_F_rows = cropped_data.F_rows/mean_F - 1
        k_kms, Pk_kms, var_kms = Pk1D.get_Pk1D(delta_F_rows,cropped_data.IVAR_rows,cropped_data.R,cropped_data.Z,z_value,z_width=0.2,N_processes=1)

        Pk1D_results += [(k_kms,Pk_kms,var_kms)]

        #Fit model
        def model_Pk_kms(k_kms,A_F,B_F):
            return tuning.P1D_z_kms_PD2013(k_kms,z_value,A_F=A_F,B_F=B_F)
        fit = curve_fit(model_Pk_kms,k_kms,Pk_kms,p0=(0.064,3.55))
        #something to do with comparing to the default values

        model_Pk_kms_fit = model_Pk_kms(k_kms,fit[0][0],fit[0][1])
        model_Pk_kms_default = model_Pk_kms(k_kms,0.064,3.55)

        """
        plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
        plt.errorbar(k_kms,Pk_kms,yerr=np.sqrt(var_kms),fmt='o',label='measured',color='orange')
        plt.plot(k_kms,model_Pk_kms_fit,label='model fit: A_F={:2.2f}, B_F={:2.2f}'.format(fit[0][0],fit[0][1]),color='b')
        plt.plot(k_kms,model_Pk_kms_default,label='model default: A_F={:2.2f}, B_F={:2.2f}'.format(0.064,3.55),color=(0.5,0.5,0.5))
        plt.fill_between(k_kms,0.9*model_Pk_kms_default,1.1*model_Pk_kms_default,label='model default +/- 10%',color=(0.5,0.5,0.5),alpha=0.5)
        #plt.plot(k_kms,independent.power_kms(z_value,k_kms,cell_size*general.get_dkms_dhMpc(z_value),False),label='added')
        plt.semilogy()
        plt.semilogx()
        plt.legend()
        plt.grid()
        plt.savefig('Pk1D_abs.pdf')
        plt.show()

        #Convert back to small cells
        #need a new function to merge cells back together

        #Measure correlation function
        cf_results += []
        files = "/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/test/[01]/*/picca-gaussian-RSD-renorm-16-*.fits"
        output = "test.fits.gz"
        rp_max = 80. #Mpc/h
        rt_max = 80. #Mpc/h
        npar = 20
        ntra = 20

        command = '/global/homes/j/jfarr/Programs/picca/bin/do_cf.py --in-dir "dummy" --from-image {} --out {} --rp-max {} --rt-max {} --np {} --nt {} --no-project --nside 64 --nproc 64'.format(files,output,rp_max,rt_max,npar,ntra)
        process = subprocess.run(command.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)

        #Run the fitter

        #For now, plot the results

        """

    return k_kms, Pk_kms, var_kms

results = []
for i in range(1):
    print('looking at:',tasks[i])
    k_kms, Pk_kms, var_kms = measure_pixel(tasks[i])
    results += [(k_kms, Pk_kms, var_kms)]

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
plt.savefig('Pk1D_abs_beta2.0.pdf')
plt.show()



# TODO: work out how to iterate
