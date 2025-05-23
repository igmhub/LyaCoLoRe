import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

N_side = 16

if len(sys.argv) > 1:
    pixels = list(range(int(sys.argv[1])))
else:
    pixels = list(range(10))

base = '../example_data/lya_skewers/'

for pixel in pixels:
    col_num = pixel/len(pixels)

    pixel100 = pixel//100
    filename = base + '/{}/{}/gaussian-colore-{}-{}.fits'.format(pixel100,pixel,N_side,pixel)

    h = fits.open(filename)

    RA = h[1].data['RA']
    DEC = h[1].data['DEC']

    plt.scatter(RA,DEC,s=2.0,c=[col_num,col_num,1-col_num])

#plt.xlim(0.,360.)
#plt.ylim(-90.,90.)
plt.xlabel('RA [deg]')
plt.ylabel('Dec [deg]')
plt.grid()
plt.show()
