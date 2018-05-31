import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

N_side = 16
pixels = list(range(4))

base = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/test/'
base = '/Users/jfarr/Projects/LyaCoLoRe/example_data/update_230518/'

for pixel in pixels:
    col_num = pixel/len(pixels)

    pixel100 = pixel//100
    filename = base + '/{}/{}/gaussian-colore-{}-{}.fits'.format(pixel100,pixel,N_side,pixel)

    h = fits.open(filename)

    RA = h[1].data['RA']
    DEC = h[1].data['DEC']

    plt.scatter(RA,DEC,s=10.0,c=[col_num,0,1-col_num])

#plt.xlim(0.,360.)
#plt.ylim(-90.,90.)
plt.grid()
plt.show()
