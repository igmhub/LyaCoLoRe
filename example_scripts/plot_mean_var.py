import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import sys
import process_functions as pf

def read_file(basedir,nside,pix):
    pix_100 = int(pix/100)
    dirname = basedir+'/'+str(pix_100)+'/'+str(pix)+'/'
    suffix = str(nside)+'-'+str(pix)+'.fits'
    # file with delta flux for picca
    filename = dirname+'/picca-gaussian-'+suffix
    #print('picca flux file',filename)
    h = fits.open(filename)
    delta_rows = h[0].data.T
    ivar_rows = h[1].data.T
    loglam = h[2].data
    h.close()
    lya=1215.67
    zs = (10**loglam)/lya - 1.0
    flux_rows = (delta_rows+1)*(pf.get_mean_F_model(zs))
    return zs,delta_rows,ivar_rows

def get_means(data_rows,ivar_rows):
    N_skewers=data_rows.shape[0]
    N_cells=data_rows.shape[1]
    print('in mean_var, we have',N_skewers,'skewers and',N_cells,'cells')
    # add tiny value to not have 0/0
    weights=np.sum(ivar_rows,axis=0)+1.e-10
    mean=np.sum(data_rows*ivar_rows,axis=0)/weights
    mean2=np.sum((data_rows**2)*ivar_rows,axis=0)/weights
    return weights,mean,mean2

# main folder where the processed files are
basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/process_output_G_hZ_4096_32_sr2.0_bm1_biasG18_picos_nside16/'
nside = 16
N_pixels = 100
pixels = np.sort(np.random.choice(list(range(12*nside**2)),size=N_pixels))

# will combine pixel statistics
sum_mean=None
sum_weights=None
sum_var=None

for pix in pixels:
    # read file for particular pixel
    zs,data_rows,ivar_rows=read_file(basedir,nside,pix)

    # compute statistics for pixel
    weights,mean,mean2 = get_means(data_rows,ivar_rows)

    # accumulate stats
    if sum_mean is None:
        sum_mean = mean*weights
        sum_mean2 = mean2*weights
        sum_weights = weights
    else:
        sum_mean += mean*weights
        sum_mean2 += mean2*weights
        sum_weights += weights

overall_mean = sum_mean/sum_weights
overall_mean2 = sum_mean2/sum_weights
overall_var = overall_mean2 - overall_mean**2
overall_sigma = np.sqrt(overall_var)

err = np.sqrt(overall_var/sum_weights)

plt.plot(zs,overall_mean,label='mean')
#plt.errorbar(zs,overall_mean,yerr=err,fmt='o',lw=2,label='mean')
#plt.plot(zs,overall_var,label='variance')
plt.plot(zs,overall_sigma,label='sigma')

plt.title('Gaussian stats: mean and variance')
plt.xlabel('z')

predicted_data = fits.open('input_files/tune_small_scale_fluctuations.fits')
z_predicted = predicted_data[1].data['z']
#mean_G_predicted = predicted_data[1].data['mean_G']
sigma_G_predicted = predicted_data[1].data['sigma_G']
#plt.plot(z_predicted,mean_G_predicted,label='predicted mean')
plt.plot(z_predicted,sigma_G_predicted,label='predicted sigma')

plt.legend()
plt.grid()
plt.savefig('mean_var_gauss.pdf')
plt.show()
