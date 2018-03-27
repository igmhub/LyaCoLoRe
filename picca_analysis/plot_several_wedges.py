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

#Set the default parameter values for wedgizing. These match picca's defaults
default_cf_parameters = {'nside': 8, 'sr': 2.0, 'rpmax': 200.0, 'rpmin': 0.0, 'rtmax': 200.0, 'np': 50, 'nt': 50, 'nr': 100, 'rmax': 200.0, 'zmin': 0.0, 'zmax': 4.0}
default_xcf_parameters = {'nside': 8, 'sr': 2.0, 'rpmax': 200.0, 'rpmin': -200.0, 'rtmax': 200.0, 'np': 100, 'nt': 50, 'nr': 100, 'rmax': 200.0, 'zmin': 0.0, 'zmax': 4.0}
file_parameters = {}

#Set the parameters to label plots by
plot_label_parameters = ['zmin','zmax']

#Either plot_per_file (i.e. each file gets its own plot, showing all mu bins)
#Or plot_per_bin (i.e. each mu bin has its own plot, but all files are grouped)
mode = 'plot_per_bin'

#Option to add in a scaled version of a CAMB power spectrum:
quantity = 'gaussian'
add_CAMB = False
CAMB_sr = ['10']
#CAMB_sr = ['10','10','10']
scale_CAMB = [1.0]
#scale_CAMB = [2.9488*0.4084/0.9998, 3.4069*0.3534/0.9998 , 3.8682*0.3110/0.9998]
#CAMB_location = '/Users/James/Projects/LyaCoLoRe/camb_scripts/'
CAMB_location = '/global/homes/j/jfarr/Projects/LyaCoLoRe/camb_scripts/'
CAMB_filename = 'camb_xi_{}.txt'

#Set up the bins of mu.
#ith bin is the range mubins[i]:mubins[i+1]
mubin_boundaries = [0.0,0.5,0.8,0.95,1.0]

mubins = []
for i in range(len(mubin_boundaries)-1):
    mubins += [(mubin_boundaries[i],mubin_boundaries[i+1])]
N_bins = len(mubins)

if mode == 'plot_per_bin':
    plot_functions.plot_per_bin(mubins,filenames,CAMB_sr)
elif mode == 'plot_per_file':
    plot_functions.plot_per_file(mubins,filenames,CAMB_sr)
else:
    error('Mode not recognised.')

plt.show()
