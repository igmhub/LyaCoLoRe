import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import healpy as hp
import glob

N_side = 16
locations = glob.glob('/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v5/test*/master.fits')

N_pixels = 12*N_side*N_side

for i,location in enumerate(locations):
    plt.figure()
    print(location)

    # identify output file we want to plot
    #hdulist = fits.open('../example_data/raw_colore_1000/out_srcs_s1_0.fits')
    hdulist = fits.open(location)

    # get quasar catalog
    initial_catalog = hdulist[1].data

    # some quasars have Dec=NaN, let's get rid of them for now
    mask = np.invert(np.isnan(initial_catalog['DEC']))
    catalog = initial_catalog[mask]
    print('initial size',len(initial_catalog))
    print('final size',len(catalog))

    # convert (RA,DEC) to (theta,phi) and find HEALPix pixel for each point
    theta,phi = np.radians(90-catalog['DEC']), np.radians(catalog['RA'])
    pixels = hp.ang2pix(N_side,theta,phi,nest=True)

    # plot angular positions of quasars color-coded by HEALPix pixel
    for pix in range(N_pixels):
        plt.scatter(catalog['RA'][pixels==pix],catalog['DEC'][pixels==pix],marker='.')

    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.xlim(0.,360.)
    plt.ylim(-90.,90.)
    plt.title(location[50:])
    plt.savefig('QSO_map_{}.png'.format(i))
    
plt.show()

