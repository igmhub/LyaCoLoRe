import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

from pyacolore import Pk1D, tuning, utils

lya = utils.lya_rest

pixels = list(range(128))
file_type = 'transmission'
z_values = [2.0,2.25,2.5,2.75,3.0,3.25]
colours = ['C0','C1','C2','C3','C4','C5']
z_width = 0.2
IVAR_cutoff = 1150.
basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v3/v3.0/'

master = fits.open(basedir + '/master.fits')
R_hMpc_mid = master[3].data['R']
z_mid = master[3].data['Z']
master.close()

k_kms_list = []
pk_kms_list = []

plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
for k,z_value in enumerate(z_values):

    print('loading data for z={}...'.format(z_value))
    skewers_list = []

    for pixel in pixels:
        print('data from pixel {}'.format(pixel),end='\r')
        pixel100 = pixel//100

        if file_type == 'transmission':
            filename = basedir + '/' + str(pixel100) + '/' + str(pixel) + '/transmission-16-' + str(pixel) + '.fits'
            h = fits.open(filename)
            lambdas = h[2].data
            z_qso = h[1].data['Z']
            pixel_skewers = h[3].data
            h.close()

            z = lambdas/lya - 1
            qsos = IVAR_cutoff*(1 + z_qso) > lya*(1 + z_value + z_width/2.)
            z_indices = abs(z - z_value) < z_width/2.

            pixel_skewers = pixel_skewers[:,z_indices]
            pixel_skewers = pixel_skewers[qsos,:]
            skewers_list += [pixel_skewers]

        elif file_type == 'picca-flux':
            filename = basedir + '/' + str(pixel100) + '/' + str(pixel) + '/picca-flux-16-' + str(pixel) + '.fits'
            h = fits.open(filename)
            log_lambdas = h[2].data
            z_qso = h[3].data['Z']
            pixel_skewers = h[0].data.T
            h.close()

            z = (10**log_lambdas)/lya - 1
            qsos = IVAR_cutoff*(1 + z_qso) > lya*(1 + z_value + z_width/2.)
            z_indices = abs(z - z_value) < z_width/2.

            pixel_skewers = pixel_skewers[:,z_indices]
            pixel_skewers = pixel_skewers[qsos,:]
            skewers_list += [pixel_skewers]

    print('')

    skewers = skewers_list[0]
    for pixel_skewers in skewers_list[1:]:
        skewers = np.concatenate((skewers,pixel_skewers),axis=0)


    if file_type == 'transmission':
        mean_F = np.average(skewers,axis=0)

        delta_F_skewers = np.zeros(skewers.shape)
        for j in range(mean_F.shape[0]):
            delta_F_skewers[:,j] = skewers[:,j]/mean_F[j] - 1

    elif file_type == 'picca-flux':
        delta_F_skewers = skewers

    R_hMpc_fine = np.interp(z,z_mid,R_hMpc_mid)[z_indices]
    z_fine = z[z_indices]

    #print(z_fine)#,delta_F_skewers.shape,z_width)
    print('calculating P1D...')
    k_kms, pk_kms, var_kms = Pk1D.get_Pk1D(delta_F_skewers,np.ones_like(delta_F_skewers),R_hMpc_fine,z_fine,z_value=z_value,z_width=z_width,units='km/s')

    k_kms_list += [k_kms]
    pk_kms_list += [pk_kms]

    plt.plot(k_kms,pk_kms,label='z={}'.format(z_value),c=colours[k])

    #Why is there a 0 k value?
    indices = k_kms > 0
    k_kms_model = np.logspace(np.log10(min(k_kms[indices])),np.log10(max(k_kms)),10**4)
    pk_kms_model = tuning.P1D_z_kms_PD2013(z_value,k_kms_model)
    plt.plot(k_kms_model,pk_kms_model,label='z={} DR9'.format(z_value),c=colours[k],linestyle=':')
    plt.fill_between(k_kms_model,pk_kms_model*1.1,pk_kms_model*0.9,color=colours[k],alpha=0.25)
    print('done!\n')

plt.axvline(x=0.005,c=(0.5,0.5,0.5))
plt.semilogx()
plt.semilogy()
plt.legend()
plt.grid()
plt.xlabel(r'k / kms$^{-1}$')
plt.ylabel('P1D')
plt.xlim(xmax=10**-1)
plt.ylim(ymin=10**0)
plt.savefig('P1D_v3.0.pdf')
plt.show()
