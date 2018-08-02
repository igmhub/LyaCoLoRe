import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

base_dir = '../example_data/lya_skewers/'
file_name = base_dir+'0/0/transmission-16-0.fits'
trans_file = fits.open(file_name)
trans_file.info()

wave = trans_file['WAVELENGTH'].data
F = trans_file['TRANSMISSION'].data
F_metal = trans_file['METALS'].data

# identify a high-z quasar to get a nice looking plot
zq = trans_file['METADATA'].data['Z']
highz = np.where(zq>3.0)[0][0]

plt.plot(wave,F[highz])
plt.xlabel('Wavelength [A]')
plt.title('Lya+Lyb transmitted flux fraction')
plt.show()

plt.plot(wave,F_metal[highz])
plt.xlabel('Wavelength [A]')
plt.title('Metal transmitted flux fraction')
plt.show()

