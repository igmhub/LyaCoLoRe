import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import process_functions as pf

alphas = np.linspace(0.0,2.0,11)
#alphas = np.linspace(0.0,20.0,11)
sigma_G_values = np.linspace(0.0,10.0,101)
beta = 1.65
Z_val = 3.0

h = fits.open('/Users/James/Projects/test_data/process_output_G_hZ_4096_32_sr2.0_bm1_nside16/nside_16_master.fits')
Z = h[2].data['Z']
D = h[2].data['D']
h.close()

D_val = np.interp(Z_val,Z,D)

plt.figure(1)
plt.figure(2)
plt.figure(3)

for alpha in alphas:
    mean_F_values = []
    sigma_dF_values = []
    sigma_flux_values = []
    for sigma_G in sigma_G_values:
        mean_F,sigma_dF = pf.get_flux_stats(sigma_G,alpha,beta,D_val)
        mean_F_values += [mean_F]
        sigma_dF_values += [sigma_dF]
        sigma_flux_values += [sigma_dF*mean_F]
    plt.figure(1)
    plt.plot(sigma_G_values,mean_F_values,label='alpha={:2.2f}'.format((alpha)))
    plt.figure(2)
    plt.plot(sigma_G_values,sigma_dF_values,label='alpha={:2.2f}'.format((alpha)))
    plt.figure(3)
    plt.plot(sigma_G_values,sigma_flux_values,label='alpha={:2.2f}'.format((alpha)))

plt.figure(1)
plt.legend()
plt.grid()
plt.xlabel('sigma_G')
plt.ylabel('mean_F')
plt.savefig('mean_F_by_sigma_G.pdf')

plt.figure(2)
plt.legend()
plt.grid()
plt.xlabel('sigma_G')
plt.ylabel('sigma_dF')
plt.savefig('sigma_dF_by_sigma_G.pdf')

plt.figure(3)
plt.legend()
plt.grid()
plt.xlabel('sigma_G')
plt.ylabel('sigma_flux')
plt.savefig('sigma_flux_by_sigma_G.pdf')

plt.show()
