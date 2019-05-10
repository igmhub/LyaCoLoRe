import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import h5py
import sys

from lyacolore import plot_functions

default_location = ['/global/homes/j/jfarr/Programs/picca/picca_analysis_000/picca_00000/']

if len(sys.argv) > 1:
    locations = sys.argv[1:]
else:
    locations = default_location

#Plotting decisions
plot_system = 'plot_per_file' #per bin or per file
#Gaussian
fit_data = {'b1': 1.0, 'b2': 1.0, 'beta1': 0., 'beta2': 0.}
#Lognormal
#fit_data = {'b1': 1.0, 'b2': 1.0, 'beta1': 0., 'beta2': 0.}
#Tau noRSD
#fit_data = {'b1': 1.65*1.13, 'b2': 1.65*1.13, 'beta1': 0., 'beta2': 0.}
#Tau
#fit_data = {'b1': 1.65*1.13, 'b2': 1.65*1.13, 'beta1': 0.9625/(1.65*1.13), 'beta2': 0.9625/(1.65*1.13)}
#fit_data = {'b1': 2.22, 'b2': 2.22, 'beta1': 0.56, 'beta2': 0.56}
#Flux
#fit_data = {'b1': -0.135, 'b2': -0.135, 'beta1': 1.178, 'beta2': 1.178}

np_bins = 40
bin_list = [0,5,10,15,20]
show_plots = True
save_plots = True
r_power = 2.

corr_objects = plot_functions.get_correlation_objects(locations)
plot_functions.make_plot_vs_rt(corr_objects,np_bins,fit_data,bin_list=bin_list,r_power=r_power,show_plots=show_plots,save_plots=save_plots)
