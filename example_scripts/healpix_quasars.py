import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import healpy as hp

# identify output file we want to plot
hdulist = fits.open('../example_data/raw_colore/N1000_out_srcs_s0_15.fits')

# get quasar catalog
initial_catalog = hdulist[1].data

# some quasars have Dec=NaN, let's get rid of them for now
mask = np.invert(np.isnan(initial_catalog['DEC']))
catalog = initial_catalog[mask]
print('initial size',len(initial_catalog))
print('final size',len(catalog))

# specify Nside parameter in HEALPix maps (power of 2)
nside=2**1

# convert (RA,DEC) to (theta,phi) and find HEALPix pixel for each point
theta,phi = np.radians(90-catalog['DEC']), np.radians(catalog['RA'])
pixels = hp.ang2pix(nside,theta,phi,nest=True)

# plot angular positions of quasars color-coded by HEALPix pixel
for pix in range(12*nside*nside):
    plt.plot(catalog['RA'][pixels==pix],catalog['DEC'][pixels==pix],'o')

plt.xlabel('RA (deg)')
plt.ylabel('Dec (deg)')
plt.show()

