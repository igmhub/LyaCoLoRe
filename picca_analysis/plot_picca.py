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
mu_boundaries = [0.0,0.5,0.8,0.95,1.0]
model = 'Slosar11'
fit_type = 'manual'
fit_data = {'b1': 1., 'b2': 1., 'beta1': 0., 'beta2': 0.}
r_power = 2
nr = 40
rmax = 160. #Mpc/h
show_plots = True
save_plots = True

rmins = [40.]#[20.,40.,60.]
afixs = ['fixed']#['free','fixed']

for rmin in rmins:
    for afix in afixs:

        rmin_label = str(rmin)
        rmin_label = rmin_label[:rmin_label.find('.')]
        suffix = '_{}r_a{}'.format(rmin_label,afix)
        res_name = '/result'+suffix+'.h5'

        #Get the correlation objects
        corr_objects = plot_functions.get_correlation_objects(locations,res_name=res_name)

        #Make plots of the objects
        plot_functions.make_plots(corr_objects,mu_boundaries,plot_system,r_power,fit_type=fit_type,fit_data=fit_data,nr=nr,rmin=rmin,rmax=rmax,save_plots=save_plots,show_plots=show_plots,suffix=suffix)
