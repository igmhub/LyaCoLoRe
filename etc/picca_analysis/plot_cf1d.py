import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import plot_functions
import utils

N_bins = 50
r_min = 0.0 #Mpc/h
r_max = 200.0 #Mpc/h
bin_size = (r_max-r_min)/N_bins

master_filepath = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/process_output_G_hZ_4096_32_sr4.0_bm1_nside8/nside_8_master.fits'

filename = sys.argv[1]
h = fits.open(filename)

DA = h[2].data['DA']

#Check that this is the numnber of contributions
N_cf1d = h[2].data['NB']

LLMIN_cf1d = h[1].header['LLMIN']
LLMAX_cf1d = h[1].header['LLMAX']
DLL_cf1d = h[1].header['DLL']
LL_cf1d = np.arange(LLMIN_cf1d,LLMAX_cf1d,DLL_cf1d)
N_cells_cf1d = LL_cf1d.shape[0]

#Get LL and R from master file
master = fits.open(master_filepath)
Z_colore = master[2].data['Z']
R_colore = master[2].data['R']
lya = utils.lya_rest
LL_colore = np.log10(lya*(1+Z_colore))

R_cf1d = np.interp(LL_cf1d,LL_colore,R_colore)

separations = abs(R_cf1d - R_cf1d[:,None])
binned_separations = (separations/bin_size).astype(int)
#binned_separations = binned_separations - (binned_separations+1)*(binned_separations>=N_bins)

r_binned = np.linspace(r_min+bin_size/2,r_max-bin_size/2,N_bins)
xi_binned = np.zeros(N_bins)
for n in range(N_bins):
    bin_entries = binned_separations==n
    N_contributions = np.sum(N_cf1d*bin_entries)
    xi_binned[n] = np.sum(bin_entries*DA*N_cf1d)/N_contributions

plt.plot(r_binned,(r_binned**2)*xi_binned)
plt.axhline(y=0,color='gray',ls=':')
plt.xlabel('r [Mpc/h]')
plt.ylabel('r^2 xi(r)')
plt.grid(True, which='both')

CAMB_filepath = '/global/homes/j/jfarr/Projects/LyaCoLoRe/camb_scripts/camb_xi_4.0.txt'
scale_CAMB = 1.0
CAMB_sr = '4.0'
r_CAMB,xi_CAMB,plot_label_CAMB = plot_functions.get_CAMB_xi(CAMB_filepath,scale_CAMB,CAMB_sr)
plt.plot(r_CAMB,xi_CAMB*(r_CAMB**2),label=plot_label_CAMB)

plt.savefig('cf1d_xir2.pdf')
plt.show()

"""
#Make a plot of data over model
xi_CAMB_interp = np.interp(r_binned,r_CAMB,xi_CAMB)
plt.plot(r_binned,xi_binned/xi_CAMB_interp,label='data/model')
plt.axhline(y=0,color='gray',ls=':')
plt.ylabel('r^2 xi(r)')
plt.grid(True, which='both')
plt.xlim(0,100)
plt.ylim(-2,3)
plt.legend()
plt.savefig('data_over_model.pdf')
plt.show()
"""
