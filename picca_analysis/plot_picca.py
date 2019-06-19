import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import h5py
import sys

from lyacolore import plot_functions

################################################################################

#Housekeeping options.
fontsize = 16
plotsize = (12, 5)
dpi = 80
show_plot = True
save_plot = True
filename = 'corr_plot.pdf'

#Create a dictionary with all information about the subplots:
# TODO: Implement a way to look up the bias and betas rather than by hand?
#'/global/homes/j/jfarr/Programs/picca/picca_analysis_042/picca_00211/'
subplots = {}
subplots[(0,0)] =  {'location':         '/Users/jfarr/Downloads/picca_00215/',
                    'filename':         'xcf_exp.fits.gz',
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax': 160.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.133, 'b2': -0.133, 'beta1': 1.4, 'beta2': 1.4},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True},
                    }
subplots[(0,1)] =  {'location':         '/Users/jfarr/Downloads/picca_00215/',
                    'filename':         'xcf_exp.fits.gz',
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax': 160.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.133, 'b2': -0.133, 'beta1': 1.4, 'beta2': 1.4},
                    'format':           {'legend': False, 'xlabel': True, 'ylabel': True},
                    }

################################################################################

#Set style options everywhere.
#plt.rc('text', usetex=True)
plt.rc('font', size=fontsize)

#Make a figure to accomodate all the subplots.
N_rows = 0
N_cols = 0
for key in subplots.keys():
    i = key[0] + 1
    j = key[1] + 1
    N_rows = np.max((i,N_rows))
    N_cols = np.max((j,N_cols))
fig, axs = plt.subplots(N_rows, N_cols, figsize=plotsize, dpi=dpi, facecolor='w', edgecolor='k')
axs = np.reshape(axs,(N_rows,N_cols))

#Make the correlation objects and plot.
for key in subplots.keys():
    corr_obj = plot_functions.get_correlation_object(subplots[key])
    subplots[key]['corr_object'] = corr_obj
    plot_functions.plot_wedges(axs[key],subplots[key])

#Save and show if desired.
plt.tight_layout()
if save_plot:
    fig.savefig(filename)
if show_plot:
    plt.show()

################################################################################

"""
for rmin in rmins:
    for afix in afixs:

        rmin_label = str(rmin)
        rmin_label = rmin_label[:rmin_label.find('.')]
        suffix = '_{}r_a{}'.format(rmin_label,afix)
        if fit_type != 'manual':
            res_name = '/result'+suffix+'.h5'
            config_name = '/config_'+suffix+'.ini'
        else:
            res_name = None

        #Get the correlation objects
        corr_objects = plot_functions.get_correlation_objects(locations,res_name=res_name)

        #Make plots of the objects
        plot_functions.make_wedge_plots(corr_objects,mu_boundaries,plot_system,r_power,fit_type=fit_type,fit_data=fit_data,nr=nr,rmin=rmin,rmax=rmax,save_plots=save_plots,show_plots=show_plots,suffix=suffix)

default_location = ['/global/homes/j/jfarr/Programs/picca/picca_analysis_000/picca_00000/']

if len(sys.argv) > 1:
    locations = sys.argv[1:]
else:
    locations = default_location

#Plotting decisions
#plot_system = 'plot_per_file' #per bin or per file
#mu_boundaries = [0.0,0.5,0.8,0.95,1.0]
#model = 'Slosar11'
#fit_type = 'picca'
#Gaussian
#fit_data = {'b1': 1.0, 'b2': 1.0, 'beta1': 0., 'beta2': 0.}
#Lognormal
#fit_data = {'b1': 1.0, 'b2': 1.0, 'beta1': 0., 'beta2': 0.}
#Tau noRSD
#fit_data = {'b1': 1.65*1., 'b2': 1.65, 'beta1': 0., 'beta2': 0.}
#Tau
#fit_data = {'b1': 1.65, 'b2': 1.65, 'beta1': 1., 'beta2': 1.}
#Flux
#fit_data = {'b1': -0.133, 'b2': -0.133, 'beta1': 1.4, 'beta2': 1.4}

#r_power = 2
#nr = 40
#rmax = 160. #Mpc/h

#rmins = [40.]#[20.,40.,60.]
#afixs = ['free','fixed']
"""
