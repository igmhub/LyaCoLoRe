import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

from pyacolore import plot_functions

default_location = ['/global/homes/j/jfarr/Programs/picca/picca_analysis_000/picca_00000/']

if len(sys.argv) > 1:
    locations = sys.argv[1:]
else:
    locations = default_location

#Plotting decisions
plot_system = 'plot_per_file' #per bin or per file
model = 'Slosar11'
fit_type = 'manual'
fit_data = {'b1': 1., 'b2': 1., 'beta1': 0., 'beta2': 0.}
r_power = 2
v_max = 10**-4
show_plots = True
save_plots = True

#Stuff we don't need for now.
rmin_label = str(rmin)
rmin_label = rmin_label[:rmin_label.find('.')]
suffix = ''
res_name = None#'/result'+suffix+'.h5'

#Get the correlation objects
corr_objects = plot_functions.get_correlation_objects(locations,res_name=res_name)

#Make plots of the objects
plot_functions.make_colour_plots(corr_objects,r_power=r_power,v_max=v_max,save_plots=True,show_plots=True,suffix=suffix)
