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
default_filename = 'cf_exp.fits.gz'
filename = sys.argv[1]

#Determine the parameters of the file
data_parameters = plot_functions.get_parameters_from_filename(filename)

#Determine the model to fit against, the parameter ranges to fit over and the assessment method
model = 'Slosar11' #'Slosar11', 'no_beta'
assessment_method = 'visual' #'visual', 'chi2'
parameter_sampling_method = 'grid' #'grid', 'MC'

from numpy import linspace as ls
b_values =      {'G': ls(1.0,1.0,1), 'D': ls(1.0,1.0,1), 'F': ls(1.0,2.0,2), 'q': ls(2.95,2.95,1)}
beta_values =   {'G': ls(0.0,0.0,1), 'D': ls(0.0,0.0,1), 'F': ls(0.0,2.0,2), 'q': ls(0.0,0.0,1)}

zmin = data_parameters['zmin']
zmax = data_parameters['zmax']
z = (zmax+zmin)/2
if zmax-zmin > 0.5:
    print('WARNING: Z RANGE IS NOT NARROW')

if assessment_method == 'visual':
    correlation_model.visual_fit(filename,b_values,beta_values,model,data_parameters,z)
    fit = []
elif assessment_method == 'chi2':
    #Maybe use the fitter2 method? Don't really understand this atm
    print('chi2 method not yet developed')
