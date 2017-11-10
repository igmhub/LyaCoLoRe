import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import healpy as hp

# identify output file we want to plot
hdulist = fits.open('../Sims/NewLya512/skewers_512_srcs_s0_0.fits')

# use only a fraction of the catalog
Nq = 100000
catalog = hdulist[1].data[0:Nq]

# specify Nside parameter in HEALPix maps (power of 2)
Nside=2**1

# convert (RA,DEC) to (theta,phi) and find HEALPix pixel for each point
pix = hp.ang2pix(Nside,(catalog['DEC']+90.0)/180.0*np.pi,catalog['RA']/180.0*np.pi)

# plot angular positions of quasars color-coded by HEALPix pixel
for i in range(12*Nside*Nside):
    plt.plot(catalog['RA'][pix==i],catalog['DEC'][pix==i],'o')

plt.xlabel('RA (deg)')
plt.ylabel('Dec (deg)')
plt.show()

