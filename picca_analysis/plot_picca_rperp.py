import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import h5py
import sys

from pyacolore import plot_functions

default_location = ['/global/homes/j/jfarr/Programs/picca/picca_analysis_000/picca_00000/']

if len(sys.argv) > 1:
    locations = sys.argv[1:]
else:
    locations = default_location

#Plotting decisions
plot_system = 'plot_per_file' #per bin or per file
fit_data = {'b1': 1., 'b2': 1., 'beta1': 0., 'beta2': 0.}
np_bins = 8
bin_list = None
show_plots = True
save_plots = True

corr_objects = plot_functions.get_correlation_objects(locations)
make_plot_vs_rt(corr_objects,np_bins,fit_data,bin_list=bin_list)
