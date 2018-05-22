import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from picca import wedgize
import sys
from numpy import linalg
import mcfit
import plot_functions
import correlation_model

#Import the picca export data
default_filenames = ['cf_exp.fits.gz']

N_files = len(sys.argv) - 1
if N_files > 0:
    filenames = sys.argv[1:]
else:
    filenames = default_filenames

print('The _exp_out file(s) to plot are:')
for filename in filenames:
    print(filename)

#Determine the model to fit against, the parameter ranges to fit over and the assessment method
model = 'Slosar11' #'Slosar11', 'no_beta'
assessment_method = 'visual' #'visual', 'chi2'
parameter_sampling_method = 'grid' #'grid', 'MC'

#Set option to ignore the below ranges and use bias from file, and beta from kaiser law
compute_b_beta = True

#Set the b and beta ranges to explore.
from numpy import linspace as ls
b_values =      {'G': ls(1.0,1.0,1), 'D': ls(1.0,1.0,1), 'F': ls(1.0,2.0,2), 'q': ls(4.333,4.333,1)}
beta_values =   {'G': ls(0.0,0.0,1), 'D': ls(0.0,0.0,1), 'F': ls(0.0,2.0,2), 'q': ls(0.225,0.225,1)}
#2.1,3.356,0.286; 2.6,4.333,0.225; 3.1,5.412,0.181; 3.6,6.589,0.150

for filename in filenames:
    #Determine the parameters of the file
    data_parameters = plot_functions.get_parameters_from_filename(filename)

    #Determine an appropriate value of z for the file.
    zmin = data_parameters['zmin']
    zmax = data_parameters['zmax']
    z = (zmax+zmin)/2
    if zmax-zmin > 0.5:
        print('WARNING: Z RANGE IS NOT NARROW')

    if assessment_method == 'visual':
        correlation_model.visual_fit(filename,b_values,beta_values,model,data_parameters,z,compute_b_beta=compute_b_beta)
        fit = []
    elif assessment_method == 'chi2':
        #Maybe use the fitter2 method? Don't really understand this atm
        print('chi2 method not yet developed')
