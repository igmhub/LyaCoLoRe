import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob

basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/process_output_G_hZsmooth_4096_32_sr2.0_bm1_biasG18_picos_newNz_mpz0_nside16/'
N_pixels = 3072

nz_files = [
            ]
nz_names = ['DESI n(z)','CoLoRe input']

fi = glob.glob(basedir+'/*/*/gaussian-colore-16-*.fits')
fi = fi[:N_pixels]

RA = []
DEC = []
Z_QSO = []

for f in fi:
    print(f)
    h = fits.open(f)
    RA = np.concatenate((RA,h[1].data['RA']]))
    DEC = np.concatenate((DEC,h[1].data['DEC']))
    Z_QSO = np.concatenate((Z_QSO,h[1].data['Z_COSMO']))
    h.close()

#m = fits.open(basedir+'/master.fits')
#m_Z_QSO = m
z_bins = np.linspace(0.,4.5,4501)
z_centres = z_bins[:-1]/2. + z_bins[1:]/2.
z_width = z_bins[1]-z_bins[0]
z_hist,_ = np.histogram(Z_QSO,bins=z_bins,density=True)

plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.bar(z_centres,z_hist,width=z_width,align='center',alpha=0.5,label='CoLoRe output')
for i,f in enumerate(nz_files):
    A = np.loadtxt(f)
    norm = np.trapz(A[:,1],A[:,0])
    plt.plot(A[:,0],A[:,1]/norm,label=nz_names[i])
plt.legend()
plt.grid()
plt.savefig('nz.pdf')
plt.show()


ang_bins = [1000,1000]
ang_hist, RA_bins, DEC_bins = np.hisogram2d(RA,DEC,bins=ang_bins)

plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.imshow(ang_hist)
plt.show()
