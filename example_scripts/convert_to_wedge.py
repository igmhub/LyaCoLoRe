import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

filename = '/global/homes/j/jfarr/Programs/picca/cf_200_GAUSS_sr2.0_bm1_quantitiesGG_RR_nside64_rpmin0.0_rpmax100.0_rtmax100.0_np25_nt25_renorm.fits.gz'
N_bins = 15
R_min = 0.
R_max = 100.
mu_bins = np.array([(0.0,0.33),(0.33,0.67),(0.67,1.0)])

R_bin_size = (R_max - R_min)/N_bins
R_edges = np.linspace(R_min,R_max,N_bins+1)

h = fits.open(filename)
RP_bins = h[1].data['RP']
RT_bins = h[1].data['RT']
DA = h[2].data['DA']
WE = h[2].data['WE']
h.close()

DA_mean = np.average(DA,weights=WE,axis=0)
DA_mean_squared = np.average(DA**2,weights=WE,axis=0)
DA_var = DA_mean_squared - DA_mean**2
WE_sum = np.sum(WE,axis=0)

R_original = np.sqrt(RP_bins**2 + RT_bins**2)
mu_original = RP_bins/R_original

R_bin_numbers = np.array(list(range(N_bins)))
binned_R_numbers = np.digitize(R_original,R_edges) - 1
binned_R_values = (R_bin_numbers + 0.5) * R_bin_size

binned_data = np.zeros((mu_bins.shape[0],R_bin_numbers.shape[0]))

for i,mu_bin in enumerate(mu_bins):
    print(mu_bin)

    mu_min = mu_bin[0]
    mu_max = mu_bin[1]

    mu_relevant = (mu_original > mu_min)*(mu_original < mu_max)

    for j,R_bin in enumerate(R_bin_numbers):
        print('   ',R_bin)

        relevant = mu_relevant * (binned_R_numbers == R_bin)
        binned_data[i,j] = np.average(DA_mean[relevant],weights=WE_sum[relevant])

    plt.scatter(binned_R_values,(binned_R_values**2)*binned_data[i,:],label=str(mu_bin))

plt.legend()
plt.grid()
plt.show()
