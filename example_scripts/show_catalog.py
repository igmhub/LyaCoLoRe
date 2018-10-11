import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob

basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/output_G_hZsmooth_4096_32_sr2.0_bm1_biasG18_picos_newNz_mpz0/'
N_files = 32
z_norm = 2.0
nz_files = ['/global/homes/j/jfarr/Projects/run_CoLoRe/input_files/Nz_qso_130618_2_colore1.txt',
            '/global/homes/j/jfarr/Projects/run_CoLoRe/input_files/Nz_qso_130618_2_colore1_hZs.txt',
            ]
nz_names = ['DESI n(z)','CoLoRe input']

fi = glob.glob(basedir+'/out_srcs_s1_*.fits')
fi = fi[:N_files]
print(len(fi))

RA = []
DEC = []
Z_QSO = []

for i,f in enumerate(fi):
    print(i+1,end='\r')
    h = fits.open(f)
    RA = np.concatenate((RA,h[1].data['RA']))
    DEC = np.concatenate((DEC,h[1].data['DEC']))
    Z_QSO = np.concatenate((Z_QSO,h[1].data['Z_COSMO']))
    h.close()

#m = fits.open(basedir+'/master.fits')
#m_Z_QSO = m
z_bins = np.linspace(0.,4.5,4501)
z_centres = z_bins[:-1]/2. + z_bins[1:]/2.
z_width = z_bins[1]-z_bins[0]
z_hist,_ = np.histogram(Z_QSO,bins=z_bins,density=False)
indices = z_centres >= z_norm
scale = np.trapz(z_hist[indices],z_centres[indices])
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.bar(z_centres,z_hist,width=z_width,align='center',alpha=0.5,label='CoLoRe output')
for i,f in enumerate(nz_files):
    A = np.loadtxt(f)
    indices = A[:,0] >= z_norm
    norm = np.trapz(A[indices,1],A[indices,0])
    plt.plot(A[:,0],scale*A[:,1]/norm,label=nz_names[i])
plt.legend(fontsize=15)
plt.grid()
plt.xlabel('z',fontsize=15)
plt.savefig('nz.pdf')
plt.show()

ang_bins = (200,100)
ang = np.array(list(zip(RA,DEC)))
ang_hist, axis = np.histogramdd(ang,bins=ang_bins)
extent = [RA.min(),RA.max(),DEC.min(),DEC.max()]
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.imshow(ang_hist,interpolation='nearest',origin='lower',extent=extent,aspect='auto')
cbar = plt.colorbar()
cbar.set_label(r'$\#$')
cbar.update_ticks()
plt.xlabel('RA',fontsize=15)
plt.ylabel('DEC',fontsize=15)
plt.savefig('ra_dec.pdf')
plt.show()
