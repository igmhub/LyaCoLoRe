import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import h5py
import sys

from pyacolore import plot_functions

default_location = ['/global/homes/j/jfarr/Programs/picca/picca_00000/']

if len(sys.argv) > 0:
    locations = sys.argv[1:]
else:
    locations = default_location

#Plotting decisions
plot_system = 'plot_per_file' #per bin or per file
mu_boundaries = [0.0,0.5,0.8,0.95,1.0]
model = 'Slosar11'
include_fits = True
r_power = 2

#Get the correlation objects
corr_objects = plot_functions.get_correlation_objects(locations)

#Make plots of the objects
plot_functions.make_plots(corr_objects,mu_boundaries,plot_system,r_power,include_fits)
