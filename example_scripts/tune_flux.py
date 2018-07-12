import numpy as np
from matplotlib.pyplot import plt
from scipy.optimize import curve_fit
import subprocess

import pixelise
import Pk1D

lambda_min = 3550.0
min_cat_z = 1.8
IVAR_cutoff = 1150.0

#Get the starting values of alpha, beta and sigma_G from file
#Decide which z values we are going to tune
z_values = [2.5]
z_width = 0.2

#Open up the Gaussian colore files
original_file_location =
input_format = 'gaussian_colore'
pixels = list(range(100))

tasks = [(pixel)]

def measure_pixel():

    location = new_base_file_location + '/' + new_file_structure.format(pixel//100,pixel)
    mean_F_data = np.array(list(zip(tuning_z_values,desired_mean_F)))

    #We work from the gaussian colore files made in 'pixelise gaussian skewers'.
    gaussian_filename = new_filename_structure.format('gaussian-colore',N_side,pixel)

    #Make a pixel object from it.
    data = pixelise.simulation_data.get_gaussian_skewers_object(location+gaussian_filename,None,input_format,SIGMA_G=measured_SIGMA_G,IVAR_cutoff=IVAR_cutoff)

    #trim skewers
    pixel_object.trim_skewers(lambda_min,min_cat_z,extra_cells=1)

    #load tuning data
    data.add_small_scale_gaussian_fluctuations(cell_size,sigma_G_z_values,extra_sigma_G_values,generator,amplitude=1.0,white_noise=False,lambda_min=0.0,IVAR_cutoff=lya)

    #Convert to flux
    data.compute_physical_skewers()
    data.compute_tau_skewers(alpha=...,beta=...)
    data.add_RSDs(alpha,beta,thermal=False)
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
        k_kms, Pk_kms, var_kms = Pk1D.get_Pk1D(skewer_rows,IVAR_rows,R_hMpc,z,z_value,z_width=0.2,N_processes=1)
        Pk1D_results += [(k_kms,Pk_kms,var_kms)]

        #Fit model
        def model_Pk_kms(k_kms,A_F,B_F):
            return Pk1D.P1D_z_kms_PD2013(k_kms,z_value,A_F=0.064,B_F=3.55)
        fit = curve_fit(model_Pk,k_kms,pk_kms)
        #something to do with comparing to the default values

        plt.plot(k_kms,Pk_kms,label='measured')
        plt.plot(k_kms,model_Pk_kms(k_kms,fit[0][0],fit[0][1]),label='model')
        plt.legend()
        plt.grid()
        plt.show()

        #Convert back to small cells
        #need a new function to merge cells back together

        #Measure correlation function
        cf_results += []
        files = "/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/test/[01]/*/picca-gaussian-RSD-renorm-16-*.fits"
        output = "test.fits.gz"
        rp_max = 80. #Mpc/h
        rt_max = 80. #Mpc/h
        np = 20
        nt = 20

        command = '/global/homes/j/jfarr/Programs/picca/bin/do_cf.py --in-dir "dummy" --from-image {} --out {} --rp-max {} --rt-max {} --np {} --nt {} --no-project --nside 64 --nproc 64'.format(files,output,rp_max,rt_max,np,nt)
        process = subprocess.run(command.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)

        #Run the fitter

        #For now, plot the results


    return



# TODO: work out how to iterate
