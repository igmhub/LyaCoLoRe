import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import process_functions as pf

alphas = np.logspace(-3.0,10.0,11)
sigma_G_values = np.linspace(0.0,10.0,101)
beta = 1.65
Z_val = 3.0

h = fits.open('/Users/jfarr/Projects/test_data/process_output_G_hZ_4096_32_sr2.0_bm1_nside16/nside_16_master.fits')
Z = h[2].data['Z']
D = h[2].data['D']
h.close()

D_val = np.interp(Z_val,Z,D)

for alpha in alphas:
    mean_F_values = []
    for sigma_G in sigma_G_values:
        mean_F,sigma_F = pf.get_flux_stats(sigma_G,alpha,beta,D_val,mean_only=True)
        mean_F_values += [mean_F]
    plt.plot(sigma_G_values,mean_F_values,label='log(alpha)={:2.2f}'.format(np.log10(alpha)))

plt.legend()
plt.grid()
plt.xlabel('sigma_G')
plt.ylabel('mean_F')
plt.savefig('mean_F_by_sigma_G_log.pdf')
plt.show()
