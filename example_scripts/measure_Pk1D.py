import numpy as np

from pyacolore import Pk1D, utils

lya = utils.lya_rest

pixels = list(range(96))
z_values = [2.0,2.25,2.5,2.75,3.0,3.25]
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

        pixel_F_skewers = h[3].data[:,z_indices][qsos,:]
        F_skewers_list += [pixel_F_skewers]

        h.close()

    print(' ')

    F_skewers = F_skewers_list[0]
    for pixel_F_skewers in F_skewers_list[1:]:
        F_skewers = np.concatenate((F_skewers,pixel_F_skewers),axis=0)

    mean_F = np.average(F_skewers,axis=1)
    delta_F_skewers = F_skewers/mean_F - 1

    R_hMpc_fine = np.interp(z,z_mid,R_hMpc_mid)
    z_fine = z

    print('calculating P1D...')
    k_kms, pk_kms, var_kms = Pk1D.get_Pk1D(delta_F_skewers,np.ones_like(delta_F_skewers),R_hMpc_fine,z_fine,z_value=z_value,z_width=z_width,units='km/s')

    k_kms_list += [k_kms]
    pk_kms_list += [pk_kms]

    plt.plot(k_kms,pk_kms,label='z='.format(z_value))
    print('done!\n')

plt.semilogx()
plt.semilogy()
plt.legend()
plt.grid()
plt.show()
