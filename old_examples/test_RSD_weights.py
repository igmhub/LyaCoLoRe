import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

from lyacolore import RSD

filename = 'example_data/lya_skewers/0/0/gaussian-colore-16-0.fits'
filename = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/process_output_G_hZsmooth_4096_32_sr2.0_bm1_biasG18_picos_newNz_mpz0_seed1003_123_nside16/0/0/gaussian-colore-16-0.fits'

h = fits.open(filename)

z_qso = h[1].data['Z_COSMO']
initial_skewers = h[2].data
velocity_skewers = h[3].data
z = h[4].data['Z']
r_hMpc = h[4].data['R']

N_qso = initial_skewers.shape[0]

first_method_final_skewers = RSD.add_skewer_RSDs(initial_skewers,None,velocity_skewers,z,r_hMpc,z_qso,thermal=False,weights=None)

new_method_weights = RSD.get_weights(None,velocity_skewers,z,r_hMpc,z_qso,thermal=False)

new_method_final_skewers = RSD.add_skewer_RSDs(initial_skewers,None,velocity_skewers,z,r_hMpc,z_qso,thermal=False,weights=new_method_weights)

print(first_method_final_skewers[0,:])
print(new_method_final_skewers[0,:])

for i in range(N_qso):
    difference = first_method_final_skewers[i,:]-new_method_final_skewers[i,:]
    if np.sum(difference) != 0:
        plt.plot(z,difference,label='difference {}'.format(i))
plt.legend()
plt.grid()
plt.show()
