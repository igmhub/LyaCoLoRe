import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import time

# TODO: make min process more efficient

output_basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/output_G_hZsmooth_4096_32_sr2.0_bm1_biasG18_picos/'
process_basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/process_output_G_hZsmooth_4096_32_sr2.0_bm1_biasG18_picos_nside16/'

N_qso = 10**5
ANG_bound = 0.2

def calculate_separations(RA,DEC):

    N_qso = RA.shape[0]

    RA_diff = np.zeros(N_qso)
    DEC_diff = np.zeros(N_qso)
    ANG_diff = np.zeros(N_qso)

    for i in range(N_qso):
        if (i+1)//100 == (i+1)/100:
            print('output file {} of {}'.format(i+1,N_qso),end='\r')

        """
        times = []
        previous_time = time.time()
        times += [time.time()-previous_time]
        previous_time = time.time()
        """

        #Centre the RA values on the current QSO, and go to angular separation
        RA_centred = RA - RA[i]
        RA_separation = 180.0 - abs(180.0 - abs(RA_centred))
        DEC_centred = DEC - DEC[i]
        DEC_separation = abs(DEC_centred)

        #Add 180 to the separation of the current QSO to rule it out
        RA_separation[i] += 180.
        DEC_separation[i] += 180.

        #Work out which QSOs are within bounds. This avoids having to calculate the angular separation of all QSOs, only those that are relatively close.
        bound = ANG_bound
        within_bound = (RA_separation<bound)*(DEC_separation<bound)

        while np.sum(within_bound) <= 1:
            bound *= 2
            print('\n  -> i={}: increasing bound from {} to {}'.format(i,bound/2.,bound))
            within_bound = [j for j in range(N_qso - 1) if RA_separation[j]<bound and DEC_separation[j]<bound]

        RA_separation = RA_separation[within_bound]
        DEC_separation = DEC_separation[within_bound]

        new_N_qso = RA_separation.shape[0]

        #Calculate the minimum separation in RA and DEC
        RA_diff[i] = np.min(RA_separation)
        DEC_diff[i] = np.min(DEC_separation)

        #Calculate the angular separation, and find the minimum
        ANG_separation = (180./np.pi)*np.arccos((np.cos((np.pi/180.)*RA_separation))*(np.cos((np.pi/180.)*DEC_separation)))
        ANG_diff[i] = np.min(ANG_separation)

    return RA_diff, DEC_diff, ANG_diff

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

#Trim off excess QSOs if necessary
if N_qso:
    RA_output = RA_output[:N_qso]
    DEC_output = DEC_output[:N_qso]
N_qso_output = RA_output.shape[0]
print('\n{} QSOs in total output'.format(N_qso_output))

print('calculating separations')
RA_output_diff, DEC_output_diff, ANG_output_diff = calculate_separations(RA_output,DEC_output)

print('\n')


#Get the RA and DEC diffs from the processed files
master = fits.open(process_basedir + '/master.fits')
RA_process = master[1].data['RA']
DEC_process = master[1].data['DEC']
master.close()

#Trim off excess QSOs if necessary
if N_qso:
    RA_process = RA_process[:N_qso]
    DEC_process = DEC_process[:N_qso]
N_qso_process = RA_process.shape[0]
print('{} QSOs in total process'.format(N_qso_process))

print('calculating separations')
RA_process_diff, DEC_process_diff, ANG_process_diff = calculate_separations(RA_process,DEC_process)
print('\n')


#Make some plots
print('calculating histograms')
hist_bins = np.linspace(0.,0.01,10001)
ang_hist_bins = np.linspace(0.,1.0,10001)
td_hist_bins = np.linspace(0.,0.0001,1001)
print('working on output RA..')
Rod_hist = np.histogram(RA_output_diff,bins=hist_bins)
print('working on output DEC..')
Dod_hist = np.histogram(DEC_output_diff,bins=hist_bins)
print('working on output ANG..')
Aod_hist = np.histogram(ANG_output_diff,bins=ang_hist_bins)
print('working on process RA..')
Rpd_hist = np.histogram(RA_process_diff,bins=hist_bins)
print('working on process DEC..')
Dpd_hist = np.histogram(DEC_process_diff,bins=hist_bins)
print('working on process ANG..')
Apd_hist = np.histogram(ANG_process_diff,bins=ang_hist_bins)
print('working on output 2D..')
tdo_hist,x_edges,y_edges = np.histogram2d(RA_output_diff,DEC_output_diff,bins=td_hist_bins)
print('working on process 2D..')
tdp_hist,x_edges,y_edges = np.histogram2d(RA_process_diff,DEC_process_diff,bins=td_hist_bins)
print('making plots...')

#1. Plot with everything wide
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.bar(Rod_hist[1][:-1], Rod_hist[0], width=np.diff(Rod_hist[1]), align="edge", label='RA output', alpha=0.5, color='r')
plt.bar(Dod_hist[1][:-1], Dod_hist[0], width=np.diff(Dod_hist[1]), align="edge", label='DEC output', alpha=0.5, color='b')
#plt.bar(Aod_hist[1][:-1], Aod_hist[0], width=np.diff(Aod_hist[1]), align="edge", label='ANG output', alpha=0.5, color='yellow')
plt.bar(Rpd_hist[1][:-1], Rpd_hist[0], width=np.diff(Rpd_hist[1]), align="edge", label='RA process', alpha=0.5, color='orange')
plt.bar(Dpd_hist[1][:-1], Dpd_hist[0], width=np.diff(Dpd_hist[1]), align="edge", label='DEC process', alpha=0.5, color='purple')
#plt.bar(Apd_hist[1][:-1], Apd_hist[0], width=np.diff(Apd_hist[1]), align="edge", label='ANG process', alpha=0.5, color='g')
plt.grid()
plt.legend()
plt.savefig('min_ang_sep_RDA_op_wide.pdf')
plt.show()

#2. Plot with everything
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.bar(Rod_hist[1][:-1], Rod_hist[0], width=np.diff(Rod_hist[1]), align="edge", label='RA output', alpha=0.5, color='r')
plt.bar(Dod_hist[1][:-1], Dod_hist[0], width=np.diff(Dod_hist[1]), align="edge", label='DEC output', alpha=0.5, color='b')
#plt.bar(Aod_hist[1][:-1], Aod_hist[0], width=np.diff(Aod_hist[1]), align="edge", label='ANG output', alpha=0.5, color='yellow')
plt.bar(Rpd_hist[1][:-1], Rpd_hist[0], width=np.diff(Rpd_hist[1]), align="edge", label='RA process', alpha=0.5, color='orange')
plt.bar(Dpd_hist[1][:-1], Dpd_hist[0], width=np.diff(Dpd_hist[1]), align="edge", label='DEC process', alpha=0.5, color='purple')
#plt.bar(Apd_hist[1][:-1], Apd_hist[0], width=np.diff(Apd_hist[1]), align="edge", label='ANG process', alpha=0.5, color='g')
plt.xlim(-0.000005,0.0001)
plt.grid()
plt.legend()
plt.savefig('min_ang_sep_RDA_op.pdf')
plt.show()

#3. Plot the two RAs
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.bar(Rod_hist[1][:-1], Rod_hist[0], width=np.diff(Rod_hist[1]), align="edge", label='RA output', alpha=0.5, color='r')
plt.bar(Rpd_hist[1][:-1], Rpd_hist[0], width=np.diff(Rpd_hist[1]), align="edge", label='RA process', alpha=0.5, color='orange')
plt.xlim(-0.000005,0.0001)
plt.grid()
plt.legend()
plt.savefig('min_ang_sep_R_op.pdf')
plt.show()

#4. Plot the two DECs
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.bar(Dod_hist[1][:-1], Dod_hist[0], width=np.diff(Dod_hist[1]), align="edge", label='DEC output', alpha=0.5, color='b')
plt.bar(Dpd_hist[1][:-1], Dpd_hist[0], width=np.diff(Dpd_hist[1]), align="edge", label='DEC process', alpha=0.5, color='purple')
plt.xlim(-0.000005,0.0001)
plt.grid()
plt.legend()
plt.savefig('min_ang_sep_D_op.pdf')
plt.show()

#5. Plot the two RAs and the two DECs
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.bar(Rod_hist[1][:-1], Rod_hist[0], width=np.diff(Rod_hist[1]), align="edge", label='RA output', alpha=0.5, color='r')
plt.bar(Rpd_hist[1][:-1], Rpd_hist[0], width=np.diff(Rpd_hist[1]), align="edge", label='RA process', alpha=0.5, color='orange')
plt.bar(Dod_hist[1][:-1], Dod_hist[0], width=np.diff(Dod_hist[1]), align="edge", label='DEC output', alpha=0.5, color='b')
plt.bar(Dpd_hist[1][:-1], Dpd_hist[0], width=np.diff(Dpd_hist[1]), align="edge", label='DEC process', alpha=0.5, color='purple')
plt.xlim(-0.000005,0.0001)
plt.grid()
plt.legend()
plt.savefig('min_ang_sep_RD_op.pdf')
plt.show()

#6. Plot the two ANGs
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.bar(Aod_hist[1][:-1], Aod_hist[0], width=np.diff(Aod_hist[1]), align="edge", label='ANG output', alpha=0.5, color='yellow')
plt.bar(Apd_hist[1][:-1], Apd_hist[0], width=np.diff(Apd_hist[1]), align="edge", label='ANG process', alpha=0.5, color='g')
plt.xlim(-0.005,0.2)
plt.grid()
plt.legend()
plt.savefig('min_ang_sep_A_op.pdf')
plt.show()

#7. Plot the output data
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.bar(Rod_hist[1][:-1], Rod_hist[0], width=np.diff(Rod_hist[1]), align="edge", label='RA output', alpha=0.5, color='r')
plt.bar(Dod_hist[1][:-1], Dod_hist[0], width=np.diff(Dod_hist[1]), align="edge", label='DEC output', alpha=0.5, color='b')
#plt.bar(Aod_hist[1][:-1], Aod_hist[0], width=np.diff(Aod_hist[1]), align="edge", label='ANG output', alpha=0.5, color='yellow')
plt.xlim(-0.000005,0.0001)
plt.grid()
plt.legend()
plt.savefig('min_ang_sep_RDA_o.pdf')
plt.show()

#8. Plot the process data
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.bar(Rpd_hist[1][:-1], Rpd_hist[0], width=np.diff(Rpd_hist[1]), align="edge", label='RA process', alpha=0.5, color='orange')
plt.bar(Dpd_hist[1][:-1], Dpd_hist[0], width=np.diff(Dpd_hist[1]), align="edge", label='DEC process', alpha=0.5, color='purple')
#plt.bar(Apd_hist[1][:-1], Apd_hist[0], width=np.diff(Apd_hist[1]), align="edge", label='ANG process', alpha=0.5, color='g')
plt.xlim(-0.000005,0.0001)
plt.grid()
plt.legend()
plt.savefig('min_ang_sep_RDA_p.pdf')
plt.show()

#9. Plot a 2D histogram of RA and DEC from the output data
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.imshow(tdo_hist,origin='low',extent=[td_hist_bins[0],td_hist_bins[-1],td_hist_bins[0],td_hist_bins[-1]])
plt.title('output')
plt.xlabel('min RA separation')
plt.ylabel('min DEC separation')
plt.savefig('min_ang_sep_2d_o.pdf')
plt.show()

#10. Plot a 2D histogram of RA and DEC from the process data
plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
plt.imshow(tdp_hist,origin='low',extent=[td_hist_bins[0],td_hist_bins[-1],td_hist_bins[0],td_hist_bins[-1]])
plt.title('process')
plt.xlabel('min RA separation')
plt.ylabel('min DEC separation')
plt.savefig('min_ang_sep_2d_o.pdf')
plt.show()

print('\nDone!')
