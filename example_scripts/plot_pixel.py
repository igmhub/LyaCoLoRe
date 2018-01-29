import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

default_pixel = 0

if len(sys.argv) > 1:
    pixel = int(sys.argv[1])
else:
    pixel = default_pixel

pixel100 = pixel//100

base = '/global/cscratch1/sd/jfarr/LyaSkewers/'
N_RA_chunks = 10

filename_basis = '{}CoLoRe_{}/process_output_4096_32/{}/{}/picca-density-8-{}.fits'
skewers_filename = filename_basis.format(base,'skewers',str(pixel100),str(pixel),str(pixel))
revamp_filename = filename_basis.format(base,'revamp',str(pixel100),str(pixel),str(pixel))

skewers = fits.open(skewers_filename)
revamp = fits.open(revamp_filename)

skewers_RA = skewers[3].data['RA']
skewers_DEC = skewers[3].data['DEC']
revamp_RA = revamp[3].data['RA']
revamp_DEC = revamp[3].data['DEC']

min_RA = 10*((min(min(skewers_RA),min(revamp_RA)))//10)
max_RA = 10*((max(max(skewers_RA),max(revamp_RA)))//10 + 1)
chunk_size = (max_RA - min_RA)/10

figs = ['']*(N_RA_chunks)

for i in range(N_RA_chunks):
    chunk_min_RA = min_RA + i*chunk_size
    chunk_max_RA = min_RA + (i+1)*chunk_size

    print('Plotting objects with RA between {} and {}.'.format(chunk_min_RA,chunk_max_RA))

    sDEC = skewers_DEC[skewers_RA >= chunk_min_RA]
    sRA = skewers_RA[skewers_RA >= chunk_min_RA]
    sDEC = sDEC[sRA < chunk_max_RA]
    sRA = sRA[sRA < chunk_max_RA]

    rDEC = revamp_DEC[revamp_RA >= chunk_min_RA]
    rRA = revamp_RA[revamp_RA >= chunk_min_RA]
    rDEC = rDEC[rRA < chunk_max_RA]
    rRA = rRA[rRA < chunk_max_RA]

    fig = plt.figure()
    a = fig.add_subplot(111)
    a.scatter(sRA,sDEC)
    a.scatter(rRA,rDEC)
    figs[i] = fig

plt.show()
