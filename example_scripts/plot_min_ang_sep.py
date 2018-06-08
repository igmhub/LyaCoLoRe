import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# TODO: Account for step from 0 to 2pi

output_basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/output_G_hZsmooth_4096_32_sr2.0_bm1_biasG18_picos/'
process_basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/process_output_G_hZsmooth_4096_32_sr2.0_bm1_biasG18_picos_nside16/'

N_qso = 6*(10**6)

#Get the RA and DEC diffs from the output files
N_output_files = 32
N_qso_output = 0
RA_output = []
DEC_output = []
for i in range(N_output_files):
    print('output file {} of {}'.format(i+1,N_output_files),end='\r')
    if N_qso_output<N_qso:
        h = fits.open(output_basedir + 'out_srcs_s1_{}.fits'.format(i))
        RA_output = np.concatenate((RA_output,h[1].data['RA']))
        DEC_output = np.concatenate((DEC_output,h[1].data['DEC']))
        N_qso_output = RA_output.shape[0]
        h.close()

if N_qso:
    RA_output = RA_output[:N_qso]
    DEC_output = DEC_output[:N_qso]

N_qso_output = RA_output.shape[0]
print('\n{} QSOs in total output'.format(N_qso_output))

RA_output_sorted = np.sort(RA_output)
DEC_output_sorted = np.sort(DEC_output)

RA_output_diff = np.zeros(N_qso_output-1)
DEC_output_diff = np.zeros(N_qso_output-1)

print('calculating separations')

RA_output_diff[0] = min(abs(RA_output_sorted[0]-RA_output_sorted[1]),abs(RA_output_sorted[0]-RA_output_sorted[-1]+360))
RA_output_diff[-1] = min(abs(RA_output_sorted[-1]-RA_output_sorted[0]-360.),abs(RA_output_sorted[-1]-RA_output_sorted[-2]))
DEC_output_diff[0] = min(abs(DEC_output_sorted[0]-DEC_output_sorted[1]),abs(DEC_output_sorted[0]-DEC_output_sorted[-1]+360))
DEC_output_diff[-1] = min(abs(DEC_output_sorted[-1]-DEC_output_sorted[0]-360.),abs(DEC_output_sorted[-1]-DEC_output_sorted[-2]))

for i in range(1,N_qso_output-2):
    if i//1000 == i/1000:
        print('output file {} of {}'.format(i+1,N_qso_output-1),end='\r')
    RA_output_diff[i] = min(abs(RA_output_sorted[i]-RA_output_sorted[i+1]),abs(RA_output_sorted[i]-RA_output_sorted[i-1]))
    DEC_output_diff[i] = min(abs(DEC_output_sorted[i]-DEC_output_sorted[i+1]),abs(DEC_output_sorted[i]-DEC_output_sorted[i-1]))

#Get the RA and DEC diffs from the processed files
master = fits.open(process_basedir + '/master.fits')

RA_process = master[1].data['RA']
DEC_process = master[1].data['DEC']

master.close()

if N_qso:
    RA_process = RA_process[:N_qso]
    DEC_process = DEC_process[:N_qso]

N_qso_process = RA_process.shape[0]
print('\n{} QSOs in total process'.format(N_qso_process))

RA_process_sorted = np.sort(RA_process)
DEC_process_sorted = np.sort(DEC_process)

RA_process_diff = np.zeros(N_qso_process-1)
DEC_process_diff = np.zeros(N_qso_process-1)

RA_process_diff[0] = min(abs(RA_process_sorted[0]-RA_process_sorted[1]),abs(RA_process_sorted[0]-RA_process_sorted[-1]+360))
RA_process_diff[-1] = min(abs(RA_process_sorted[-1]-RA_process_sorted[0]-360.),abs(RA_process_sorted[-1]-RA_process_sorted[-2]))
DEC_process_diff[0] = min(abs(DEC_process_sorted[0]-DEC_process_sorted[1]),abs(DEC_process_sorted[0]-DEC_process_sorted[-1]+360))
DEC_process_diff[-1] = min(abs(DEC_process_sorted[-1]-DEC_process_sorted[0]-360.),abs(DEC_process_sorted[-1]-DEC_process_sorted[-2]))

print('calculating separations')
for i in range(1,N_qso_process-2):
    if i//1000 ==i/1000:
        print('process file {} of {}'.format(i+1,N_qso_process-1),end='\r')
    RA_process_diff[i] = min(abs(RA_process_sorted[i]-RA_process_sorted[i+1]),abs(RA_process_sorted[i]-RA_process_sorted[i-1]))
    DEC_process_diff[i] = min(abs(DEC_process_sorted[i]-DEC_process_sorted[i+1]),abs(DEC_process_sorted[i]-DEC_process_sorted[i-1]))


#Make some plots
hist_bins=np.linspace(0.,0.01,10001)

#1. Plot with everything
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
Rod_hist=plt.hist(RA_output_diff,bins=hist_bins,label='RA output',alpha=0.5)
Dod_hist=plt.hist(DEC_output_diff,bins=hist_bins,label='DEC output',alpha=0.5)
Rpd_hist=plt.hist(RA_process_diff,bins=hist_bins,label='RA process',alpha=0.5)
Dpd_hist=plt.hist(DEC_process_diff,bins=hist_bins,label='DEC process',alpha=0.5)
plt.xlim(-0.000005,0.0001)
plt.grid()
plt.legend()
plt.savefig('min_ang_sep_RD_op.pdf')
#plt.show()

#2. Plot the two RAs
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
Rod_hist=plt.hist(RA_output_diff,bins=hist_bins,label='RA output',alpha=0.5)
Rpd_hist=plt.hist(RA_process_diff,bins=hist_bins,label='RA process',alpha=0.5)
plt.xlim(-0.000005,0.0001)
plt.grid()
plt.legend()
plt.savefig('min_ang_sep_R_op.pdf')
#plt.show()

#3. Plot the two DECs
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
Dod_hist=plt.hist(DEC_output_diff,bins=hist_bins,label='DEC output',alpha=0.5)
Dpd_hist=plt.hist(DEC_process_diff,bins=hist_bins,label='DEC process',alpha=0.5)
plt.xlim(-0.000005,0.0001)
plt.grid()
plt.legend()
plt.savefig('min_ang_sep_D_op.pdf')
#plt.show()

#4. Plot the output data
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
Rod_hist=plt.hist(RA_output_diff,bins=hist_bins,label='RA output',alpha=0.5)
Dod_hist=plt.hist(DEC_output_diff,bins=hist_bins,label='DEC output',alpha=0.5)
plt.xlim(-0.000005,0.0001)
plt.grid()
plt.legend()
plt.savefig('min_ang_sep_RD_o.pdf')
#plt.show()

#5. Plot the process data
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
Rpd_hist=plt.hist(RA_process_diff,bins=hist_bins,label='RA process',alpha=0.5)
Dpd_hist=plt.hist(DEC_process_diff,bins=hist_bins,label='DEC process',alpha=0.5)
plt.xlim(-0.000005,0.0001)
plt.grid()
plt.legend()
plt.savefig('min_ang_sep_RD_p.pdf')
#plt.show()

print('\nDone!')
