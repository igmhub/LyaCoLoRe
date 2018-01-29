import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import lya_mock_functions as mock
import sys

base_dir ='/Users/jfarr/Projects/repixelise/test_output/test_multi/'
h = fits.open(base_dir + 'nside_8_master.fits')
pixel_list = list(sorted(set(h[1].data['PIXNUM'])))

for pixel in pixel_list:
    pixel_100 = int(pixel/100)
    picca_density = fits.open(base_dir + str(pixel_100) + '/' + str(pixel) + '/' + 'picca-density-8-' + str(pixel) + '.fits')
    picca_flux = fits.open(base_dir + str(pixel_100) + '/' + str(pixel) + '/' + 'picca-flux-8-' + str(pixel) + '.fits')

    sum_IVAR_density = sum(IVAR_density)
    sum_IVAR_flux = sum(IVAR_flux)

    if sum_IVAR_density == 0:
        print('Pixel {}:'.format(pixel))
        print(' -> dens skewer ivar min = {}'.format(min(sum_IVAR_density)))
    elif sum_IVAR_flux == 0:
        print('Pixel {}:'.format(pixel))
        print(' -> flux skewer ivar min = {}'.format(min(sum_IVAR_flux)))
    else:
        print('Pixel {} ivar ok.'.format(pixel))
