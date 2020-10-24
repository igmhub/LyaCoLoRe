import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

from lyacolore import plot_functions

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
vmax = 10**-4
show_plots = True
save_plots = True

filename = 'berkeley_cross_0.2_DLA.pdf'
figsize=(12,8)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_045/picca_00240/',
                    'filename':         'xcf_exp_0.2_randoms.fits.gz',
                    #'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1049, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': False, 'xlabel': True, 'ylabel': True},
                    }


#Get the correlation objects
corr_object = plot_functions.get_correlation_object(subplots[(0,0)])

#Make plots of the objects
plot_functions.make_colour_plots(corr_object,vmax=vmax,save_plots=True,show_plots=True,suffix=suffix)
