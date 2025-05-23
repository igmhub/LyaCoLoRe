import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import h5py
import sys

from lyacolore import plot_functions

################################################################################

#Housekeeping options.
fontsize = 18
figsize = (12, 5)
dpi = 80
show_plot = True
save_plot = True
filename = 'corr_plot_systematics.pdf'

#Create a dictionary with all information about the subplots:

#Main method correlations plot:
filename = 'berkeley_cross_0.2_DLA_vs_rt.pdf'
figsize=(12,8)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_045/picca_00240/',
                    'filename':         'xcf_exp_0.2_randoms.fits.gz',
                    'rp_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0),(12.0,16.0)],
                    'rp_bin_colours':   ['C0','C1','C2','C3'],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True},
                    }

"""
filename = 'corr_plot_vs_rt.pdf'
subplots = {}
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_037/combined/',
                    'filename':         'cf_exp.fits.gz',
                    'rp_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0),(60.0,64.0)],
                    'rp_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax': 160.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.133, 'b2': -0.133, 'beta1': 1.4, 'beta2': 1.4},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True},
                    }
subplots[(0,1)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_041/picca_00223/',
                    'filename':         'xcf_exp_2400k_shuffle.fits.gz',
                    'rp_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0),(60.0,64.0)],
                    'rp_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax': 160.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.133, 'b2': 2.0, 'beta1': 1.4, 'beta2': 0.79},
                    'format':           {'legend': False, 'xlabel': True, 'ylabel': True},
                    }

"""
"""

#Systematics correlations plot:
filename = 'corr_plot_systematics_vs_rt.pdf'
subplots = {}
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_042/picca_00224/',
                    'filename':         'cf_exp_800k.fits.gz',
                    'rp_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0),(60.0,64.0)],
                    'rp_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax': 160.0},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1049, 'b2': -0.1049, 'beta1': 1.3783, 'beta2': 1.3783},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True},
                    }
subplots[(0,1)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_041/picca_00214/',
                    'filename':         'xcf_exp_noshuffle.fits.gz',
                    'rp_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0),(60.0,64.0)],
                    'rp_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax': 160.0},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1049, 'b2': 2.0, 'beta1': 1.4, 'beta2': 0.79},
                    'format':           {'legend': False, 'xlabel': True, 'ylabel': True},
                    }
"""
"""

#HCDs test plot:
filename = 'corr_plot_HCDs_vs_rt.pdf'
subplots = {}
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_041/picca_00226/',
                    'filename':         'xcf_exp_2400k_pdcov.fits.gz',
                    'rp_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0),(60.0,64.0)],
                    'rp_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax': 160.0},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1049, 'b2': 2.0, 'beta1': 1.3783, 'beta2': 0.79},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True},
                    }
subplots[(0,1)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_041/picca_00227/',
                    'filename':         'xcf_exp_zb0.05_2400k_pdcov.fits.gz',
                    'rp_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0),(60.0,64.0)],
                    'rp_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax': 160.0},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1049, 'b2': 2.0, 'beta1': 1.3783, 'beta2': 0.79},
                    'format':           {'legend': False, 'xlabel': True, 'ylabel': True},
                    }
"""

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
fig, axs = plt.subplots(N_rows, N_cols, figsize=figsize, dpi=dpi, facecolor='w', edgecolor='k')
axs = np.reshape(axs,(N_rows,N_cols))

#Make the correlation objects and plot.
for key in subplots.keys():
    corr_obj = plot_functions.get_correlation_object(subplots[key])
    subplots[key]['corr_object'] = corr_obj
    plot_functions.plot_rp_bins_vs_rt(axs[key],subplots[key])

#Save and show if desired.
plt.tight_layout()
if save_plot:
    fig.savefig(filename)
if show_plot:
    plt.show()

################################################################################
