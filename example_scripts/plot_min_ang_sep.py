import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

output_basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/output_G_hZsmooth_4096_32_sr2.0_bm1_biasG18_picos'
process_basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/process_output_G_hZsmooth_4096_32_sr2.0_bm1_biasG18_picos_nside16'

#Get the RA and DEC diffs from the output files
N_output_files = 32
RA_output = []
DEC_output = []
for i in range(N_output_files):
    h = fits.open(output_basedir + 'out_srcs_s1_{}.fit'.format(i))
    RA_output = np.concatenate((RA_output,h[1].data['RA']))
    DEC_output = np.concatenate((DEC_output,h[1].data['DEC']))

RA_output_sorted = np.sort(RA_output)
DEC_output_sorted = np.sort(DEC_output)

RA_output_diff = np.zeros(N_qso-1)
DEC_output_diff = np.zeros(N_qso-1)

for i in range(N_qso-1):
    RA_output_diff[i] = abs(RA_output_sorted[i]-RA_output_sorted[i+1])
    DEC_output_diff[i] = abs(DEC_output_sorted[i]-DEC_output_sorted[i+1])


#Get the RA and DEC diffs from the processed files
master = fits.open(process_basedir + '/master.fits')

RA_process = master[1].data['RA']
DEC_process = master[1].data['DEC']
N_qso = RA.shape[0]

RA_process_sorted = np.sort(RA_process)
DEC_process_sorted = np.sort(DEC_process)

RA_process_diff = np.zeros(N_qso-1)
DEC_process_diff = np.zeros(N_qso-1)

for i in range(N_qso-1):
    RA_process_diff[i] = abs(RA_process_sorted[i]-RA_process_sorted[i+1])
    DEC_process_diff[i] = abs(DEC_process_sorted[i]-DEC_process_sorted[i+1])


#Make some plots
hist_bins=np.linspace(0.,0.01,10001)

#1. Plot with everything
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
Rod_hist=plt.hist(RA_output_diff,bins=hist_bins,label='RA output',alpha=0.5)
Dod_hist=plt.hist(DEC_output_diff,bins=hist_bins,label='DEC output',alpha=0.5)
Rpd_hist=plt.hist(RA_process_diff,bins=hist_bins,label='RA process',alpha=0.5)
Dpd_hist=plt.hist(DEC_process_diff,bins=hist_bins,label='DEC process',alpha=0.5)
plt.xlim(-0.00001,0.0001)
plt.grid()
plt.legend()
plt.show()

#2. Plot the two RAs
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
Rod_hist=plt.hist(RA_output_diff,bins=hist_bins,label='RA output',alpha=0.5)
Rpd_hist=plt.hist(RA_process_diff,bins=hist_bins,label='RA process',alpha=0.5)
plt.xlim(-0.00001,0.0001)
plt.grid()
plt.legend()
plt.show()

#3. Plot the two DECs
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
Dod_hist=plt.hist(DEC_output_diff,bins=hist_bins,label='DEC output',alpha=0.5)
Dpd_hist=plt.hist(DEC_process_diff,bins=hist_bins,label='DEC process',alpha=0.5)
plt.xlim(-0.00001,0.0001)
plt.grid()
plt.legend()
plt.show()

#4. Plot the output data
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
Rod_hist=plt.hist(RA_output_diff,bins=hist_bins,label='RA output',alpha=0.5)
Dod_hist=plt.hist(DEC_output_diff,bins=hist_bins,label='DEC output',alpha=0.5)
plt.xlim(-0.00001,0.0001)
plt.grid()
plt.legend()
plt.show()

#5. Plot the process data
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
Rpd_hist=plt.hist(RA_process_diff,bins=hist_bins,label='RA process',alpha=0.5)
Dpd_hist=plt.hist(DEC_process_diff,bins=hist_bins,label='DEC process',alpha=0.5)
plt.xlim(-0.00001,0.0001)
plt.grid()
plt.legend()
plt.show()
