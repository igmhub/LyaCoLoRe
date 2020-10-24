import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

if len(sys.argv) > 1:
    location = sys.argv[1]
else:
    location = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/process_output_G_hZ_4096_32_sr2.0_bm1_nside8/'
h = fits.open(location+'nside_8_master.fits')

try:
    Z_QSO_NO_RSD = h[1].data['Z_QSO_NO_RSD']
    Z_QSO_RSD = h[1].data['Z_QSO_RSD']
    plt.hist(Z_QSO_NO_RSD,bins=500,label='Z_QSO_NO_RSD')
    #plt.hist(Z_QSO_RSD,bins=500,label='Z_QSO_RSD')
except KeyError:
    Z_QSO = h[1].data['Z_QSO']
    plt.hist(Z_QSO,bins=500,label='Z_QSO')

plt.xlim(1.8,3.9)
plt.legend()
plt.grid()

plt.savefig('catalog_histogram.pdf')
plt.show()
