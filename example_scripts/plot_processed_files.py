import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import sys

# main folder where the processed files are
basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_revamp/process_output_hZ_4096_32/'
nside = 8

pix = 3
#pix = int(sys.argv[1])
pix_100 = int(pix/100)

dirname = basedir+'/'+str(pix_100)+'/'+str(pix)+'/'
suffix = str(nside)+'-'+str(pix)+'.fits'
print('dir name',dirname)

# file with physical density
phys_dens_file = dirname+'/physical-colore-'+suffix
print('physical density file',phys_dens_file)
phys_dens_hdu = fits.open(phys_dens_file)
phys_dens_hdu.info()
zq = phys_dens_hdu[1].data['Z_COSMO']
phys_dens = phys_dens_hdu[2].data
phys_dens_zs = phys_dens_hdu[4].data['Z']
phys_dens_hdu.close()

# file with Gaussian density
gauss_dens_file = dirname+'/gaussian-colore-'+suffix
print('gaussian density file',gauss_dens_file)
gauss_dens_hdu = fits.open(gauss_dens_file)
gauss_dens_hdu.info()
gauss_dens = gauss_dens_hdu[2].data
gauss_dens_zs = gauss_dens_hdu[4].data['Z']
gauss_dens_hdu.close()

# file with delta gaussian for picca
picca_gauss_file = dirname+'/picca-gaussian-'+suffix
print('picca gaussian file',picca_gauss_file)
picca_gauss_hdu = fits.open(picca_gauss_file)
picca_gauss_hdu.info()
picca_gauss_delta = picca_gauss_hdu[0].data
picca_ivar = picca_gauss_hdu[1].data
picca_gauss_ll = picca_gauss_hdu[2].data
picca_gauss_hdu.close()
lya=1215.67
picca_zs = (10**picca_gauss_ll)/lya - 1.0

# file with delta density for picca
picca_dens_file = dirname+'/picca-density-'+suffix
print('picca density file',picca_dens_file)
picca_dens_hdu = fits.open(picca_dens_file)
picca_dens_hdu.info()
picca_dens_delta = picca_dens_hdu[0].data
picca_dens_hdu.close()

# plot all fields, for a given skewer
iskewer=0
plt.plot(phys_dens_zs,phys_dens[iskewer],lw=2,label='physical density')
plt.plot(gauss_dens_zs,gauss_dens[iskewer],lw=2,label='gauss density')
plt.plot(picca_zs,picca_ivar[:,iskewer],label='ivar picca')
plt.plot(picca_zs,picca_gauss_delta[:,iskewer],ls=':',label='gauss picca')
plt.plot(picca_zs,picca_dens_delta[:,iskewer],ls=':',label='density picca')
plt.axvline(x=zq[iskewer])
plt.legend()
plt.show()

