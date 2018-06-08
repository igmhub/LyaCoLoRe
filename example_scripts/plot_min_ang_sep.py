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
ANG_output_diff = np.zeros(N_qso_output-1)

print('calculating separations')

"""
RA_output_diff[0] = min(abs(RA_output_sorted[0]-RA_output_sorted[1]),abs(RA_output_sorted[0]-RA_output_sorted[-1]+360))
RA_output_diff[-1] = min(abs(RA_output_sorted[-1]-RA_output_sorted[0]-360.),abs(RA_output_sorted[-1]-RA_output_sorted[-2]))
DEC_output_diff[0] = min(abs(DEC_output_sorted[0]-DEC_output_sorted[1]),abs(DEC_output_sorted[0]-DEC_output_sorted[-1]+360))
DEC_output_diff[-1] = min(abs(DEC_output_sorted[-1]-DEC_output_sorted[0]-360.),abs(DEC_output_sorted[-1]-DEC_output_sorted[-2]))

for i in range(1,N_qso_output-2):
    if i//1000 == i/1000:
        print('output file {} of {}'.format(i+1,N_qso_output-1),end='\r')
    RA_output_diff[i] = min(abs(RA_output_sorted[i]-RA_output_sorted[i+1]),abs(RA_output_sorted[i]-RA_output_sorted[i-1]))
    DEC_output_diff[i] = min(abs(DEC_output_sorted[i]-DEC_output_sorted[i+1]),abs(DEC_output_sorted[i]-DEC_output_sorted[i-1]))
print('\n')
"""

for i in range(N_qso_output - 1):
    if (i+1)//1000 == (i+1)/1000:
        print('output file {} of {}'.format(i+1,N_qso_output-1),end='\r')
    RA_output_centred = RA_output - RA_output[i]
    RA_output_separation = 180.0 - abs(180.0 - abs(RA_output_centred))
    DEC_output_centred = RA_output - RA_output[i]
    DEC_output_separation = abs(RA_output_centred)

    RA_output_diff[i] = np.min(RA_output_separation)
    DEC_output_diff[i] = np.min(DEC_output_separation)

    ANG_output_separation = np.arccos((cos(RA_output_separation))*(cos(DEC_output_separation)))
    ANG_output_diff[i] = np.min(ANG_output_separation)
print('\n')


#Get the RA and DEC diffs from the processed files
master = fits.open(process_basedir + '/master.fits')

RA_process = master[1].data['RA']
DEC_process = master[1].data['DEC']

master.close()

if N_qso:
    RA_process = RA_process[:N_qso]
    DEC_process = DEC_process[:N_qso]

N_qso_process = RA_process.shape[0]
print('{} QSOs in total process'.format(N_qso_process))

RA_process_sorted = np.sort(RA_process)
DEC_process_sorted = np.sort(DEC_process)

RA_process_diff = np.zeros(N_qso_process-1)
DEC_process_diff = np.zeros(N_qso_process-1)
ANG_process_diff = np.zeros(N_qso_process-1)

"""
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
print('\n')
"""

for i in range(N_qso_process - 1):
    if (i+1)//1000 == (i+1)/1000:
        print('process file {} of {}'.format(i+1,N_qso_process-1),end='\r')
    RA_process_centred = RA_process - RA_process[i]
    RA_process_separation = 180.0 - abs(180.0 - abs(RA_process_centred))
    DEC_process_centred = RA_process - RA_process[i]
    DEC_process_separation = abs(RA_process_centred)

    RA_process_diff[i] = np.min(RA_process_separation)
    DEC_process_diff[i] = np.min(DEC_process_separation)

    ANG_process_separation = np.arccos((cos(RA_process_separation))*(cos(DEC_process_separation)))
    ANG_process_diff[i] = np.min(ANG_process_separation)
print('\n')


#Make some plots
print('calculating histograms')
hist_bins=np.linspace(0.,0.01,10001)
print('working on output RA..')
Rod_hist=np.histogram(RA_output_diff,bins=hist_bins)
print('working on output DEC..')
Dod_hist=np.histogram(DEC_output_diff,bins=hist_bins)
print('working on output ANG..')
Aod_hist=np.histogram(ANG_output_diff,bins=hist_bins)
print('working on process RA..')
Rpd_hist=np.histogram(RA_process_diff,bins=hist_bins)
print('working on process DEC..')
Dpd_hist=np.histogram(DEC_process_diff,bins=hist_bins)
print('working on process ANG..')
Apd_hist=np.histogram(ANG_process_diff,bins=hist_bins)
print('making plots...')

#1. Plot with everything
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.bar(Rod_hist[1][:-1], Rod_hist[0], width=np.diff(Rod_hist[1]), align="edge", label='RA output', alpha=0.5)
plt.bar(Dod_hist[1][:-1], Dod_hist[0], width=np.diff(Dod_hist[1]), align="edge", label='DEC output', alpha=0.5)
plt.bar(Aod_hist[1][:-1], Aod_hist[0], width=np.diff(Aod_hist[1]), align="edge", label='ANG output', alpha=0.5)
plt.bar(Rpd_hist[1][:-1], Rpd_hist[0], width=np.diff(Rpd_hist[1]), align="edge", label='RA process', alpha=0.5)
plt.bar(Dpd_hist[1][:-1], Dpd_hist[0], width=np.diff(Dpd_hist[1]), align="edge", label='DEC process', alpha=0.5)
plt.bar(Apd_hist[1][:-1], Apd_hist[0], width=np.diff(Apd_hist[1]), align="edge", label='ANG process', alpha=0.5)
plt.xlim(-0.000005,0.0001)
plt.grid()
plt.legend()
plt.savefig('min_ang_sep_RDA_op.pdf')
plt.show()

#2. Plot the two RAs
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.bar(Rod_hist[1][:-1], Rod_hist[0], width=np.diff(Rod_hist[1]), align="edge", label='RA output', alpha=0.5)
plt.bar(Rpd_hist[1][:-1], Rpd_hist[0], width=np.diff(Rpd_hist[1]), align="edge", label='RA process', alpha=0.5)
plt.xlim(-0.000005,0.0001)
plt.grid()
plt.legend()
plt.savefig('min_ang_sep_R_op.pdf')
plt.show()

#3. Plot the two DECs
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.bar(Dod_hist[1][:-1], Dod_hist[0], width=np.diff(Dod_hist[1]), align="edge", label='DEC output', alpha=0.5)
plt.bar(Dpd_hist[1][:-1], Dpd_hist[0], width=np.diff(Dpd_hist[1]), align="edge", label='DEC process', alpha=0.5)
plt.xlim(-0.000005,0.0001)
plt.grid()
plt.legend()
plt.savefig('min_ang_sep_D_op.pdf')
plt.show()

#4. Plot the two ANGs
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.bar(Aod_hist[1][:-1], Aod_hist[0], width=np.diff(Aod_hist[1]), align="edge", label='ANG output', alpha=0.5)
plt.bar(Apd_hist[1][:-1], Apd_hist[0], width=np.diff(Apd_hist[1]), align="edge", label='ANG process', alpha=0.5)
plt.xlim(-0.000005,0.0001)
plt.grid()
plt.legend()
plt.savefig('min_ang_sep_A_op.pdf')
plt.show()

#5. Plot the output data
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.bar(Rod_hist[1][:-1], Rod_hist[0], width=np.diff(Rod_hist[1]), align="edge", label='RA output', alpha=0.5)
plt.bar(Dod_hist[1][:-1], Dod_hist[0], width=np.diff(Dod_hist[1]), align="edge", label='DEC output', alpha=0.5)
plt.bar(Aod_hist[1][:-1], Aod_hist[0], width=np.diff(Aod_hist[1]), align="edge", label='ANG output', alpha=0.5)
plt.xlim(-0.000005,0.0001)
plt.grid()
plt.legend()
plt.savefig('min_ang_sep_RDA_o.pdf')
plt.show()

#6. Plot the process data
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.bar(Rpd_hist[1][:-1], Rpd_hist[0], width=np.diff(Rpd_hist[1]), align="edge", label='RA process', alpha=0.5)
plt.bar(Dpd_hist[1][:-1], Dpd_hist[0], width=np.diff(Dpd_hist[1]), align="edge", label='DEC process', alpha=0.5)
plt.bar(Apd_hist[1][:-1], Apd_hist[0], width=np.diff(Apd_hist[1]), align="edge", label='ANG process', alpha=0.5)
plt.xlim(-0.000005,0.0001)
plt.grid()
plt.legend()
plt.savefig('min_ang_sep_RDA_p.pdf')
plt.show()

print('\nDone!')
