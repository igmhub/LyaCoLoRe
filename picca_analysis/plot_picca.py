import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import h5py
import sys

from lyacolore import plot_functions

################################################################################

#Housekeeping options.
fontsize = 18
dpi = 80
show_plot = True
save_plot = True
basedir = '/project/projectdirs/desi/users/jfarr/LyaCoLoRe_paper/'

#Create a dictionary with all information about the subplots:

"""
figsize = (12,6)
filename = 'compare_lya_qso_cross_covariances_stack.pdf'
subplots = {}
subplots[(0,0)] =  {'location':         basedir+'/analysis/correlation_functions/stack/measurements/lya_qso_cross/',
                    'filename':         '/correlations/xcf_exp_lya_qso_cross.fits.gz',
                    #'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_qso_cross/',
                    #'result_name':      'result_lya_qso_cross_rmin40.0_rmax160.0_afree_smoothingfree.h5',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_auto__lya_qso_cross/',
                    'result_name':      'result_lya_auto__lya_qso_cross_rmin40.0_rmax160.0_afree.h5',
                    'corr_name':        'LYA(LYA)xQSO',
                    'abs_mu':           True,
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.133, 'b2': 2.0, 'beta1': 1.4, 'beta2': 0.79},
                    'format':           {'legend': False, 'xlabel': True, 'ylabel': True,'title': 'SS+r', 'leg_loc': 'shared'},
                    }
subplots[(0,1)] =  {'location':         basedir+'/analysis/correlation_functions/stack/measurements/lya_qso_cross/',
                    'filename':         '/correlations/xcf_exp_lya_qso_cross_subsamples_nosmooth.fits.gz',
                    #'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_qso_cross/',
                    #'result_name':      'result_lya_qso_cross_rmin40.0_rmax160.0_afree_smoothingfree.h5',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_auto__lya_qso_cross/',
                    'result_name':      'result_lya_auto__lya_qso_cross_rmin40.0_rmax160.0_afree_subsamples_nosmooth.h5',
                    'corr_name':        'LYA(LYA)xQSO',
                    'abs_mu':           True,
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.133, 'b2': 2.0, 'beta1': 1.4, 'beta2': 0.79},
                    'format':           {'legend': False, 'xlabel': True, 'ylabel': True,'title': 'All SS, no smooth', 'leg_loc': 'shared'},
                    }
"""
"""
figsize = (18,6)
filename = 'compare_lya_qso_cross_covariances_v9.0.0.pdf'
subplots = {}
subplots[(0,0)] =  {'location':         basedir+'/analysis/correlation_functions/v9.0.0/measurements/lya_qso_cross/',
                    'filename':         '/correlations/xcf_exp_lya_qso_cross.fits.gz',
                    #'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_qso_cross/',
                    #'result_name':      'result_lya_qso_cross_rmin40.0_rmax160.0_afree_smoothingfree.h5',
                    'result_location':  basedir+'/analysis/correlation_functions/v9.0.0/fits/lya_auto__lya_qso_cross/',
                    'result_name':      'result_lya_auto__lya_qso_cross_rmin40.0_rmax160.0_afree.h5',
                    'corr_name':        'LYA(LYA)xQSO',
                    'abs_mu':           True,
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.133, 'b2': 2.0, 'beta1': 1.4, 'beta2': 0.79},
                    'format':           {'legend': False, 'xlabel': True, 'ylabel': True,'title': 'SS+r', 'leg_loc': 'shared'},
                    }
subplots[(0,1)] =  {'location':         basedir+'/analysis/correlation_functions/v9.0.0/measurements/lya_qso_cross/',
                    'filename':         '/correlations/xcf_exp_lya_qso_cross_pdcov.fits.gz',
                    #'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_qso_cross/',
                    #'result_name':      'result_lya_qso_cross_rmin40.0_rmax160.0_afree_smoothingfree.h5',
                    'result_location':  basedir+'/analysis/correlation_functions/v9.0.0/fits/lya_auto__lya_qso_cross/',
                    'result_name':      'result_lya_auto__lya_qso_cross_rmin40.0_rmax160.0_afree_pdcov.h5',
                    'corr_name':        'LYA(LYA)xQSO',
                    'abs_mu':           True,
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.133, 'b2': 2.0, 'beta1': 1.4, 'beta2': 0.79},
                    'format':           {'legend': False, 'xlabel': True, 'ylabel': True,'title': 'pdcov', 'leg_loc': 'shared'},
                    }
subplots[(0,2)] =  {'location':         basedir+'/analysis/correlation_functions/v9.0.0/measurements/lya_qso_cross/',
                    'filename':         '/correlations/xcf_exp_lya_qso_cross_wick_0.99_10000.fits.gz',
                    #'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_qso_cross/',
                    #'result_name':      'result_lya_qso_cross_rmin40.0_rmax160.0_afree_smoothingfree.h5',
                    'result_location':  basedir+'/analysis/correlation_functions/v9.0.0/fits/lya_auto__lya_qso_cross/',
                    'result_name':      'result_lya_auto__lya_qso_cross_rmin40.0_rmax160.0_afree_wick_0.99_10000.h5',
                    'corr_name':        'LYA(LYA)xQSO',
                    'abs_mu':           True,
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.133, 'b2': 2.0, 'beta1': 1.4, 'beta2': 0.79},
                    'format':           {'legend': False, 'xlabel': True, 'ylabel': True,'title': 'Wick', 'leg_loc': 'shared'},
                    }
"""
"""
figsize = (12,6)
filename = 'lya_dla_cross_rmin40_joint_mupm_subbin.pdf'
subplots = {}
subplots[(0,0)] =  {'location':         basedir+'/analysis/correlation_functions/stack/measurements/lya_dla_cross/',
                    'filename':         '/correlations/xcf_exp_lya_dla_cross_subbin.fits.gz',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_auto__lya_dla_cross/',
                    'result_name':      'result_lya_auto__lya_dla_cross_rmin40.0_rmax160.0_afree_subbin.h5',
                    'corr_name':        'LYA(LYA)xDLA',
                    'abs_mu':           False,
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 50, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.121, 'b2': 2.17, 'beta1': 1.6, 'beta2': 0.446},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'leg_loc': 0},
                    }
subplots[(0,1)] =  {'location':         basedir+'/analysis/correlation_functions/stack/measurements/lya_dla_cross/',
                    'filename':         '/correlations/xcf_exp_lya_dla_cross_subbin.fits.gz',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_auto__lya_dla_cross/',
                    'result_name':      'result_lya_auto__lya_dla_cross_rmin40.0_rmax160.0_afree_subbin.h5',
                    'corr_name':        'LYA(LYA)xDLA',
                    'abs_mu':           False,
                    'mu_bins':          [(-0.5,0.0),(-0.8,-0.5),(-0.95,-0.8),(-1.0,-0.95)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 50, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.121, 'b2': 2.17, 'beta1': 1.6, 'beta2': 0.446},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'leg_loc': 0},
                    }
"""
"""
#Main method correlations plot.
figsize = (12,6)
filename = 'corr_plot_rmin40_joint.pdf'
subplots = {}
subplots[(0,0)] =  {'location':         basedir+'/analysis/correlation_functions/stack/measurements/lya_auto/',
                    'filename':         '/correlations/cf_exp_lya_auto.fits.gz',
                    #'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_auto/',
                    #'result_name':      'result_lya_auto_rmin40.0_rmax160.0_afree_smoothingfree.h5',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_auto__lya_qso_cross/',
                    'result_name':      'result_lya_auto__lya_qso_cross_rmin40.0_rmax160.0_afree.h5',
                    'corr_name':        'LYA(LYA)xLYA(LYA)',
                    'abs_mu':           True,
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.1049, 'b2': -0.1049, 'beta1': 1.3783, 'beta2': 1.3783},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'leg_loc': 'shared'},
                    }
subplots[(0,1)] =  {'location':         basedir+'/analysis/correlation_functions/stack/measurements/lya_qso_cross/',
                    'filename':         '/correlations/xcf_exp_lya_qso_cross.fits.gz',
                    #'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_qso_cross/',
                    #'result_name':      'result_lya_qso_cross_rmin40.0_rmax160.0_afree_smoothingfree.h5',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_auto__lya_qso_cross/',
                    'result_name':      'result_lya_auto__lya_qso_cross_rmin40.0_rmax160.0_afree.h5',
                    'corr_name':        'LYA(LYA)xQSO',
                    'abs_mu':           True,
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.133, 'b2': 2.0, 'beta1': 1.4, 'beta2': 0.79},
                    'format':           {'legend': False, 'xlabel': True, 'ylabel': True, 'leg_loc': 'shared'},
                    }
"""
"""
figsize = (12,6)
filename = 'corr_plot_systematics_rmin20_metalsv9.0.0.pdf'
subplots = {}
subplots[(0,0)] =  {'location':         basedir+'/analysis/correlation_functions/v9.0.0/measurements/lya_aa_auto/',
                    'filename':         '/correlations/cf_exp_lya_aa_auto.fits.gz',
                    'result_location':  basedir+'/analysis/correlation_functions/v9.0.0/fits/lya_aa_auto/',
                    'result_name':      'result_lya_aa_auto_rmin20.0_rmax160.0_afree.h5',
                    'corr_name':        'LYA(LYA)xLYA(LYA)',
                    'abs_mu':           True,
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 20., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.1049, 'b2': -0.1049, 'beta1': 1.3783, 'beta2': 1.3783},
                    'fit_plot_data':    {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'leg_loc': 'shared'},
                    }
subplots[(0,1)] =  {'location':         basedir+'/analysis/correlation_functions/stack/measurements/lya_dla_cross/',
                    'filename':         '/correlations/xcf_exp_lya_dla_cross.fits.gz',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_dla_cross/',
                    'result_name':      'result_lya_dla_cross_rmin40.0_rmax160.0_afree.h5',
                    'corr_name':        'LYA(LYA)xDLA',
                    'abs_mu':           True,
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.133, 'b2': 2.0, 'beta1': 1.4, 'beta2': 0.79},
                    'fit_plot_data':    {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'format':           {'legend': False, 'xlabel': True, 'ylabel': True, 'leg_loc': 'shared'},
                    }
"""
figsize = (12,6)
filename = 'corr_plot_systematics_rmin40_sep_temp_dla_auto.pdf'
subplots = {}
subplots[(0,0)] =  {'location':         basedir+'/analysis/correlation_functions/v9.0.0/measurements/lya_aa_auto/',
                    'filename':         '/correlations/cf_exp_lya_aa_auto.fits.gz',
                    'result_location':  basedir+'/analysis/correlation_functions/v9.0.0/fits/lya_aa_auto/',
                    'result_name':      'result_lya_aa_auto_rmin20.0_rmax160.0_afree.h5',
                    'corr_name':        'LYA(LYA)xLYA(LYA)',
                    'abs_mu':           True,
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 20., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.1049, 'b2': -0.1049, 'beta1': 1.3783, 'beta2': 1.3783},
                    'fit_plot_data':    {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'leg_loc': 'shared'},
                    }
subplots[(0,1)] =  {'location':         basedir+'/analysis/correlation_functions/stack/measurements/dla_auto/',
                    'filename':         '/correlations/co_exp_dla_auto.fits.gz',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/dla_auto/',
                    'result_name':      'result_dla_auto_rmin40.0_rmax160.0_afree.h5',
                    'corr_name':        'DLAxDLA',
                    'abs_mu':           True,
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,1.0)],
                    #'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 80.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.133, 'b2': 2.0, 'beta1': 1.4, 'beta2': 0.79},
                    'format':           {'legend': False, 'xlabel': True, 'ylabel': True, 'leg_loc': 'shared'},
                    }
"""
figsize = (12,6)
filename = 'corr_plot_systematics_rmin40_sep_temp_dla_auto.pdf'
subplots = {}
subplots[(0,0)] =  {'location':         basedir+'/analysis/correlation_functions/stack/measurements/dla_auto/',
                    'filename':         '/correlations/co_exp_dla_auto.fits.gz',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/dla_auto/',
                    'result_name':      'result_dla_auto_rmin40.0_rmax160.0_afree.h5',
                    'corr_name':        'DLAxDLA',
                    'abs_mu':           True,
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,1.0)],
                    #'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 50, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.133, 'b2': 2.0, 'beta1': 1.4, 'beta2': 0.79},
                    'format':           {'legend': False, 'xlabel': True, 'ylabel': True, 'title': 'f fixed', 'leg_loc': 'shared'},
                    }
subplots[(0,1)] =  {'location':         basedir+'/analysis/correlation_functions/stack/measurements/dla_auto/',
                    'filename':         '/correlations/co_exp_dla_auto.fits.gz',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/dla_auto/',
                    'result_name':      'result_dla_auto_rmin40.0_rmax160.0_afree_ffree.h5',
                    'corr_name':        'DLAxDLA',
                    'abs_mu':           True,
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,1.0)],
                    #'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 50, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.133, 'b2': 2.0, 'beta1': 1.4, 'beta2': 0.79},
                    'format':           {'legend': False, 'xlabel': True, 'ylabel': True,'title': 'f free', 'leg_loc': 'shared'},
                    }
"""
"""
figsize = (12,6)
filename = 'corr_plot_qso_auto_rmin40.pdf'
subplots = {}
subplots[(0,0)] =  {'location':         basedir+'/analysis/correlation_functions/stack/measurements/qso_auto/',
                    'filename':         'correlations/co_exp_qso_auto.fits.gz',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/qso_auto/',
                    'result_name':      'result_qso_auto_rmin40.0_rmax160.0_afree.h5',
                    'corr_name':        'QSOxQSO',
                    'abs_mu':           True,
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.1049, 'b2': -0.1049, 'beta1': 1.3783, 'beta2': 1.3783},
                    'fit_plot_data':    {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'format':           {'legend': True, 'title': None, 'xlabel': True, 'ylabel': True, 'leg_loc': 'shared'},
                    }
"""
"""
figsize = (12,6)
filename = 'lya_dla_compare_methods.pdf'
subplots = {}
subplots[(0,0)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_061/picca_00322/',
                    'filename':         'xcf_exp_randoms.fits.gz',
                    'abs_mu':           True,
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 50, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 20., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.1049, 'b2': -0.1049, 'beta1': 1.3783, 'beta2': 1.3783},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'leg_loc': 'shared', 'title': r'$\lambda_{r,DLA}<1100 \AA $'},
                    }
subplots[(0,1)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_061/picca_00322/',
                    'filename':         'xcf_exp_randoms.fits.gz',
                    'abs_mu':           True,
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 50, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.133, 'b2': 2.0, 'beta1': 1.4, 'beta2': 0.79},
                    'format':           {'legend': False, 'xlabel': True, 'ylabel': True, 'leg_loc': 'shared', 'title': r'$\mu<0,\ \lambda_{r,DLA}<1200 \AA$'},
                    }
"""
"""
figsize=(12,8)
filename = 'dla_auto_0.5.pdf'
subplots = {}
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_046/picca_00251/',
                    'filename':         'co_exp.fits.gz',
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,1.0)],
                    #'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 20., 'rmax': 160., 'afix': 'fixed'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': 2.0, 'b2': 2.0, 'beta1': 0.48, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': None, 'leg_loc': 0},
                    }
"""
"""
filename = 'berkeley_auto_0.2.pdf'
subplots = {}
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_043/picca_00229/',
                    'filename':         'cf_exp_0.2.fits.gz',
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1246, 'b2': -0.1246, 'beta1': 1.5343, 'beta2': 1.5343},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True},
                    }
"""
"""
filename = 'corr_plot_metals.pdf'
figsize = (12,8)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_045/picca_00239/',
                    'filename':         'cf_exp_0.2.fits.gz',
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 160.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 10., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.1049, 'b2': -0.1049, 'beta1': 1.3783, 'beta2': 1.3783},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': None},
                    }
"""
"""
filename = 'lya_qso_cross_0.2.pdf'
figsize=(15,6)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_043/picca_00230/',
                    'filename':         'xcf_exp_0.2_shuffle.fits.gz',
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'abs_mu':           True,
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 20., 'rmax': 160., 'afix': 'fixed'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': 'Lya x QSO'},
                    }
subplots[(0,1)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_043/picca_00230/',
                    'filename':         'xcf_exp_0.2_shuffle.fits.gz',
                    'mu_bins':          [(-0.5,0.0),(-0.8,-0.5),(-0.95,-0.8),(-1.0,-0.95)],
                    'abs_mu':           True,
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 20., 'rmax': 160., 'afix': 'fixed'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': 'Lya x QSO'},
                    }

"""
"""
filename = 'lya_dla_cross_0.2_lr1150.pdf'
figsize=(15,6)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_047/picca_00252/',
                    'filename':         'xcf_exp_0.2_randoms.fits.gz',
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'abs_mu':           False,
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'fixed'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': 'Lya x DLA', 'leg_loc': 0},
                    }
subplots[(0,1)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_047/picca_00252/',
                    'filename':         'xcf_exp_0.2_randoms.fits.gz',
                    'mu_bins':          [(-0.5,0.0),(-0.8,-0.5),(-0.95,-0.8),(-1.0,-0.95)],
                    'abs_mu':           False,
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'fixed'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': 'Lya x DLA', 'leg_loc': 0},
                    }
"""
"""
filename = 'lya_dla_cross_0.2_subbin.pdf'
figsize=(15,6)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_045/picca_00248/',
                    'filename':         'xcf_exp_0.2_randoms_subbin.fits.gz',
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'abs_mu':           False,
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'fixed'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': '', 'leg_loc': 0},
                    }
subplots[(0,1)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_045/picca_00248/',
                    'filename':         'xcf_exp_0.2_randoms_subbin.fits.gz',
                    'mu_bins':          [(-0.5,0.0),(-0.8,-0.5),(-0.95,-0.8),(-1.0,-0.95)],
                    'abs_mu':           False,
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'fixed'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': '', 'leg_loc': 0},
                    }
"""
"""
filename = 'berkeley_cross_0.2_DLA.pdf'
figsize=(12,8)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_045/picca_00240/',
                    'filename':         'xcf_exp_0.2_randoms.fits.gz',
                    #'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bins':          [(0.95,1.0),(-1.0,-0.95)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 200.0},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True},
                    }

"""
"""
#Systematics correlations plot:
figsize=(12,6)
filename = 'corr_plot_systematics.pdf'
subplots = {}
subplots[(0,0)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_045/picca_00239/',
                    'filename':         'cf_exp_0.2.fits.gz',
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'abs_mu':           False,
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 160.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 20., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.1049, 'b2': -0.1049, 'beta1': 1.3783, 'beta2': 1.3783},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'leg_loc': 'shared'},
                    }
subplots[(0,1)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_045/picca_00240/',
                    'filename':         'xcf_exp_0.2_randoms_subbin.fits.gz',
                    'abs_mu':           True,
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 160.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 20., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.1049, 'b2': 2.0, 'beta1': 1.4, 'beta2': 0.79},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': '', 'leg_loc': 'shared'},
                    }
"""
"""
#HCDs test plot:
filename = 'corr_plot_HCDs.pdf'
subplots = {}
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_045/picca_00246/',
                    'filename':         'xcf_exp_0.2_randoms.fits.gz',
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 160.0},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1049, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True},
                    }
subplots[(0,1)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_045/picca_00240/',
                    'filename':         'xcf_exp_0.2_randoms.fits.gz',
                    'mu_bins':          [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.0)],
                    'mu_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax_plot': 160.0},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax':160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1049, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
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
share_legend = False
for key in subplots.keys():
    corr_obj = plot_functions.get_correlation_object(subplots[key])
    if subplots[key]['format']['leg_loc'] == 'shared':
        share_legend = True
    subplots[key]['corr_object'] = corr_obj
    plot_functions.plot_wedges(fig,axs[key],subplots[key])

#Save and show if desired.
if share_legend:
    rect = (0,0.1,1,1.0)
else:
    rect=(0,0.0,1,1.0)
plt.tight_layout(rect=rect)
plt.margins(x=0.1, y=0.1)
if save_plot:
    fig.savefig(filename)
if show_plot:
    plt.show()

################################################################################
