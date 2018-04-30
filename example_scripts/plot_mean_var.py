import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import sys

def read_file(basedir,nside,pix):
    pix_100 = int(pix/100)
    dirname = basedir+'/'+str(pix_100)+'/'+str(pix)+'/'
    suffix = str(nside)+'-'+str(pix)+'.fits'
    # file with delta flux for picca
    filename = dirname+'/picca-flux-'+suffix
    print('picca flux file',filename)
    h = fits.open(filename)
    delta_rows = h[0].data.T
    ivar_rows = h[1].data.T
    loglam = h[2].data
    h.close()
    lya=1215.67
    zs = (10**loglam)/lya - 1.0
    return zs,delta_rows,ivar_rows

def get_means(delta_rows,ivar_rows):
    N_skewers=delta.shape[0]
    N_cells=delta.shape[1]
    print('in mean_var, we have',N_skewers,'skewers and',N_cells,'cells')
    # add tiny value to not have 0/0
    weights=np.sum(ivar,axis=0)+1.e-10
    mean=np.sum(delta*ivar,axis=0)/weights
    mean2=np.sum((delta**2)*ivar,axis=1)/weights
    return weights,mean,mean2

# main folder where the processed files are
basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_revamp/process_output_hZ_4096_32/'
nside = 16
N_pixels = 10
pixels = np.sort(np.random.choice(list(range(12*nside**2)),size=N_pixels))

# will combine pixel statistics
sum_mean=None
sum_weights=None
sum_var=None

for pix in pixels:
    # read file for particular pixel
    zs,delta_rows,ivar_rows=read_file(basedir,nside,pix)

    # compute statistics for pixel
    weights,mean,mean2 = get_mean_var(delta,ivar)

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

err = np.sqrt(sum_var/sum_weights)

plt.errorbar(zs,overall_mean,yerr=err,fmt='o',lw=2,label='mean')
plt.plot(zs,overall_var,lw=2,label='variance')

plt.title('Flux delta stats: mean and variance')
plt.xlabel('z')

plt.legend()
plt.grid()
plt.savefig('mean_var_gauss.pdf')
plt.show()
