import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import healpy as hp

# identify output file we want to plot
hdulist = fits.open('../example_data/raw_colore/test_N1000.fits')

# get quasar catalog
initial_catalog = hdulist[1].data

# some quasars have Dec=NaN, let's get rid of them for now
mask = np.invert(np.isnan(initial_catalog['DEC']))
catalog = initial_catalog[mask]
print('initial size',len(initial_catalog))
print('final size',len(catalog))

# specify Nside parameter in HEALPix maps (power of 2)
Nside=2**2

# convert (RA,DEC) to (theta,phi) and find HEALPix pixel for each point
pix = hp.ang2pix(Nside,(catalog['DEC']+90.0)/180.0*np.pi,catalog['RA']/180.0*np.pi)

# plot angular positions of quasars color-coded by HEALPix pixel
for i in range(12*Nside*Nside):
    plt.plot(catalog['RA'][pix==i],catalog['DEC'][pix==i],'o')

plt.xlabel('RA (deg)')
plt.ylabel('Dec (deg)')
plt.show()

