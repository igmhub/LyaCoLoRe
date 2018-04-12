import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

N_bins = 50
r_min = 0.0 #Mpc/h
r_max = 200.0 #Mpc/h
bin_size = (r_max-r_min)/N_bins

master_filepath = ''

filename = sys.argv[1]
h = fits.open(filename)

DA = h[2].data['DA']

#Check that this is the numnber of contributions
N_cf1d = h[2].data['NB']

LLMIN_cf1d = h[1].header['LLMIN']
LLMAX_cf1d = h[1].header['LLMAX']
DLL_cf1d = h[1].header['DLL']
LL_cf1d = np.arange(LLMIN,LLMAX,DLL)
N_cells_cf1d = LL_cf1d.shape[0]

#Get LL and R from master file
master = fits.open(master_filepath)
Z_colore = master[2].data['Z']
R_colore = master[2].data['R']
lya = 1215.67
LL_colore = np.log10(lya*(1+Z))

R_cf1d = np.interp(LL_cf1d,LL_colore,R_colore)

separations = abs(R_cf1d - R_cf1d[:,None])
binned_separations = (separations/bin_size).astype(int)
#binned_separations = binned_separations - (binned_separations+1)*(binned_separations>=N_bins)

R_binned = np.linspace(r_min+bin_size/2,r_max-bin_size/2,N_bins)
xi_binned = np.zeros(N_bins)
for n in range(N_bins):
    bin_entries = binned_separations==n
    N_contributions = np.sum(N_cf1d*bin_entries)
    xi_binned[n] = np.sum(bin_entries*DA*N_cf1d)/N_contributions

plt.plot(R_binned,xi_binned)
plt.show()
