import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import sys

def read_file(basedir,nside,pix):
    pix_100 = int(pix/100)
    dirname = basedir+'/'+str(pix_100)+'/'+str(pix)+'/'
    suffix = str(nside)+'-'+str(pix)+'.fits'
    # file with delta gaussian for picca
    filename = dirname+'/picca-gaussian-'+suffix
    print('picca gaussian file',filename)
    h = fits.open(filename)
    delta = h[0].data
    ivar = h[1].data
    loglam = h[2].data
    h.close()
    lya=1215.67
    zs = (10**loglam)/lya - 1.0
    return zs,delta,ivar

def get_mean_var(delta,ivar):
    Ns=delta.shape[0]
    Np=delta.shape[1]
    print('in mean_var, we have',Ns,'skewers and',Np,'pixels')
    # add tiny value to not have 0/0
    weights=np.sum(ivar,axis=1)+1.e-10
    mean=np.sum(delta*ivar,axis=1)/weights
    var=np.sum((delta**2)*ivar,axis=1)/weights
    return weights,mean,var

# main folder where the processed files are
basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_revamp/process_output_hZ_4096_32/'
nside = 8
pixels = [0,100,200,300,400,500,600,700]

# first figure out Gaussian smoothing
filename=basedir+'0/0/physical-colore-'+str(nside)+'-0.fits'
h = fits.open(filename)
h.info()
h[4].header
sigma_g=h[4].header['SIGMA_G']
print('sigma=',sigma_g)

# will combine pixel statistics
sum_mean=None
sum_weights=None
sum_var=None

for pix in pixels:
    # read file for particular pixel
    zs,delta,ivar=read_file(basedir,nside,pix)

    # compute statistics for pixel
    weights,mean,var = get_mean_var(delta,ivar)
    err_mean=np.sqrt(var/weights)
    plt.errorbar(zs,mean,yerr=err_mean,fmt='o')
    plt.plot(zs,var)

    # accumulate stats
    if sum_mean is None:
        sum_mean = mean*weights
        sum_var = var*weights
        sum_weights = weights
    else:
        sum_mean += mean*weights
        sum_var += var*weights
        sum_weights += weights

sum_mean /= sum_weights
sum_var /= sum_weights
sum_err = np.sqrt(sum_var/sum_weights)
plt.errorbar(zs,sum_mean,yerr=sum_err,fmt='o',lw=2,color='black')
plt.plot(zs,sum_var,lw=2,color='black')
plt.ylim(-0.2,1.5)
plt.axhline(y=0,lw=0.5,ls='--',color='gray')
plt.axhline(y=sigma_g**2,lw=0.5,ls='--',color='gray')
plt.title('Gaussian delta stats: mean and variance')
plt.xlabel('z')
plt.savefig('mean_var_gauss.pdf')
plt.show()

