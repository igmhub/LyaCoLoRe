import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from picca import wedgize
import sys
from numpy import linalg
import plot_functions

default_filenames = ['cf_exp.fits.gz']

#Set up the list of files to wedgize.
N_files = len(sys.argv) - 1
if N_files > 0:
    filenames = sys.argv[1:]
else:
    filenames = default_filenames

print('The cf_exp_out file(s) to wedgize are:')
for filename in filenames:
    print(filename)

#Set the parameters to label plots by
plot_label_parameters = ['correlation','quantity','sr']

#Either plot_per_file (i.e. each file gets its own plot, showing all mu bins)
#Or plot_per_bin (i.e. each mu bin has its own plot, but all files are grouped)
mode = 'plot_per_bin'

#Option to add in a scaled version of a CAMB power spectrum:
quantity = 'gaussian'
add_CAMB = True
#CAMB_sr = ['10']
CAMB_sr = ['10','20','40']
scale_CAMB = [1.0,1.0,1.0]
#scale_CAMB = [2.9488*0.4084/0.9998, 3.4069*0.3534/0.9998 , 3.8682*0.3110/0.9998]
#CAMB_location = '/Users/James/Projects/LyaCoLoRe/camb_scripts/'
CAMB_location = '/global/homes/j/jfarr/Projects/LyaCoLoRe/camb_scripts/'
CAMB_filename = 'camb_xi_{}.txt'

#Set up the bins of mu.
#ith bin is the range mubins[i]:mubins[i+1]
mubin_boundaries = [-1.0,1.0]

mubins = []
for i in range(len(mubin_boundaries)-1):
    mubins += [(mubin_boundaries[i],mubin_boundaries[i+1])]
N_bins = len(mubins)

if mode == 'plot_per_bin':
    plot_functions.plot_per_bin(mubins,filenames,add_CAMB,plot_label_parameters,CAMB_sr=CAMB_sr,scale_CAMB=scale_CAMB)
elif mode == 'plot_per_file':
    plot_functions.plot_per_file(mubins,filenames,add_CAMB,plot_label_parameters,CAMB_sr=CAMB_sr,scale_CAMB=scale_CAMB)
else:
    error('Mode not recognised.')

plt.show()
