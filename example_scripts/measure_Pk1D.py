import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

from pyacolore import Pk1D, tuning, utils

lya = utils.lya_rest

pixels = list(range(10))
z_values = [2.0,2.5,3.0]
z_width = 0.2
IVAR_cutoff = 1150.
basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/process_output_G_hZsmooth_4096_32_sr2.0_bm1_biasG18_picos_newNz_mpz0_nside16/'

master = fits.open(basedir + '/master.fits')
R_hMpc_mid = master[3].data['R']
z_mid = master[3].data['Z']
master.close()

k_kms_list = []
pk_kms_list = []

plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
for z_value in z_values:

    print('loading data for z={}...'.format(z_value))
    F_skewers_list = []

    for pixel in pixels:
        print('data from pixel {}'.format(pixel),end='\r')
        pixel100 = pixel//100
        filename = basedir + '/' + str(pixel100) + '/' + str(pixel) + '/transmission-16-' + str(pixel) + '.fits'
        h = fits.open(filename)
        lambdas = h[2].data
        z_qso = h[1].data['Z_noRSD']

        z = lambdas/lya - 1
        qsos = IVAR_cutoff*(1 + z_qso) > lya*(1 + z_value + z_width/2.)
        z_indices = abs(z - z_value) < z_width/2.

        pixel_F_skewers = h[3].data
        pixel_F_skewers = pixel_F_skewers[:,z_indices]
        pixel_F_skewers = pixel_F_skewers[qsos,:]
        F_skewers_list += [pixel_F_skewers]

        h.close()

    print(' ')

    F_skewers = F_skewers_list[0]
    for pixel_F_skewers in F_skewers_list[1:]:
        F_skewers = np.concatenate((F_skewers,pixel_F_skewers),axis=0)

    mean_F = np.average(F_skewers,axis=0)

    delta_F_skewers = np.zeros(F_skewers.shape)
    for j in range(mean_F.shape[0]):
        delta_F_skewers[:,j] = F_skewers[:,j]/mean_F[j] - 1

    R_hMpc_fine = np.interp(z,z_mid,R_hMpc_mid)[z_indices]
    z_fine = z[z_indices]

    #print(z_fine)#,delta_F_skewers.shape,z_width)
    print('calculating P1D...')
    k_kms, pk_kms, var_kms = Pk1D.get_Pk1D(delta_F_skewers,np.ones_like(delta_F_skewers),R_hMpc_fine,z_fine,z_value=z_value,z_width=z_width,units='km/s')

    k_kms_list += [k_kms]
    pk_kms_list += [pk_kms]

    plt.plot(k_kms,pk_kms,label='z={}'.format(z_value))
    
    #Why is there a 0 k value?
    indices = k_kms > 0
    k_kms_model = np.logspace(np.log10(min(k_kms[indices])),np.log10(max(k_kms)),10**4)
    pk_kms_model = tuning.P1D_z_kms_PD2013(z_value,k_kms_model)
    plt.plot(k_kms_model,pk_kms_model,label='z={} DR9'.format(z_value))
    print('done!\n')

plt.semilogx()
plt.semilogy()
plt.legend()
plt.grid()
plt.show()
