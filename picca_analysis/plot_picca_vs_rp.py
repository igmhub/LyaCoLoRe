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
basedir = '/project/projectdirs/desi/users/jfarr/LyaCoLoRe_paper/'

"""
filename = 'lya_qso_cross_vs_rp_singlebin.pdf'
figsize = (12,6)
subplots = {}
subplots[(0,0)] =  {'location':         basedir+'/analysis/correlation_functions/stack/measurements/lya_qso_cross/',
                    'filename':         '/correlations/xcf_exp_lya_qso_cross.fits.gz',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_auto__lya_qso_cross/',
                    'result_name':      'result_lya_auto__lya_qso_cross_rmin40.0_rmax160.0_afree.h5',
                    'corr_name':        'LYA(LYA)xQSO',
                    'rt_bins':          [(0.0,4.0)],
                    #'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0)],
                    'rt_bin_colours':   ['k','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.6, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': None, 'leg_loc': 'shared'},
                    }
"""
filename = 'corr_plot_vs_rp_lya_hcd_cross_stack_cbf.pdf'
figsize = (12,6)
subplots = {}
subplots[(0,0)] =  {'location':         basedir+'/analysis/correlation_functions/stack/measurements/lya_dla_cross/',
                    'filename':         '/correlations/xcf_exp_lya_dla_cross_norandoms.fits.gz',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_auto__lya_dla_cross/',
                    'result_name':      'result_lya_auto__lya_dla_cross_rmin40.0_rmax160.0_afree.h5',
                    'corr_name':        'LYA(LYA)xDLA',
                    #'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0)],
                    #'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'rt_bin_colours':   ['#F5793A','#A95AA1','#85C0F9','#0F2080'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.6, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$\lambda_{r,\rm{HCD}}<1100\ \rm{\AA}$', 'leg_loc': 'shared'},
                    }
subplots[(0,1)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_053/picca_00281/',
#subplots[(0,1)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_054/picca_00286/',
                    'filename':         'xcf_exp_randoms.fits.gz',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_auto__lya_dla_cross/',
                    'result_name':      'result_lya_auto__lya_dla_cross_rmin40.0_rmax160.0_afree.h5',
                    'corr_name':        'LYA(LYA)xDLA',
                    #'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0)],
                    #'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'rt_bin_colours':   ['#F5793A','#A95AA1','#85C0F9','#0F2080'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.6, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$\lambda_{r,\rm{HCD}}<\lambda_{\rm{Ly}\alpha}$', 'leg_loc': 'shared'},
#                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$\lambda_{r,\rm{HCD}}<1200\ \rm{\AA}$', 'leg_loc': 'shared'},
                    }
"""
filename = 'lya_dla_cross_shuffled_forests.pdf'
figsize=(16,8)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_062/picca_00325/',
                    'filename':         '/xcf_exp_norandoms.fits.gz',
                    'corr_name':        'LYA(LYA)xDLA',
                    #'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.6, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': 'data only', 'leg_loc': 'shared'},
                    }
subplots[(0,1)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_062/picca_00325/',
                    'filename':         '/xcf_exp_randomsonly.fits.gz',
                    'corr_name':        'LYA(LYA)xDLA',
                    #'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.6, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': 'randoms only', 'leg_loc': 'shared'},
}
"""
"""
filename = 'lya_dla_cross_testnzqso_contributions.pdf'
figsize=(16,16)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_062/picca_00323/',
                    'filename':         '/xcf_exp_norandoms.fits.gz',
                    'corr_name':        'LYA(LYA)xDLA',
                    #'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.6, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': '59 Q/sqd, data only', 'leg_loc': 'shared'},
                    }
subplots[(1,0)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_062/picca_00323/',
                    'filename':         '/xcf_exp_randomsonly.fits.gz',
                    'corr_name':        'LYA(LYA)xDLA',
                    #'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.6, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': '59 Q/sqd, randoms only', 'leg_loc': 'shared'},
                    }
subplots[(0,1)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_062/picca_00324/',
                    'filename':         '/xcf_exp_norandoms.fits.gz',
                    'corr_name':        'LYA(LYA)xDLA',
                    #'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.6, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': '118 Q/sqd, data only', 'leg_loc': 'shared'},
                    }
subplots[(1,1)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_062/picca_00324/',
                    'filename':         '/xcf_exp_randomsonly.fits.gz',
                    'corr_name':        'LYA(LYA)xDLA',
                    #'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.6, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': '118 Q/sqd, randoms only', 'leg_loc': 'shared'},
                    }

"""
"""
filename = 'lya_dla_cross_randomsonly_lrcuts_altbins.pdf'
figsize=(16,16)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_053/picca_00281/',
                    'filename':         'xcf_exp_randomsonly.fits.gz',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_auto__lya_dla_cross/',
                    'result_name':      'result_lya_auto__lya_dla_cross_rmin40.0_rmax160.0_afree.h5',
                    'corr_name':        'LYA(LYA)xDLA',
                    #'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.6, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$\lambda_{r,\rm{DLA}}<\lambda_{\rm{Ly}\alpha}$', 'leg_loc': 'shared'},
                    }
subplots[(0,1)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_054/picca_00286/',
                    'filename':         'xcf_exp_randomsonly.fits.gz',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_auto__lya_dla_cross/',
                    'result_name':      'result_lya_auto__lya_dla_cross_rmin40.0_rmax160.0_afree.h5',
                    'corr_name':        'LYA(LYA)xDLA',
                    #'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.6, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$\lambda_{r,\rm{DLA}}<1200 \AA$', 'leg_loc': 'shared'},
                    }
subplots[(1,0)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_055/picca_00291/',
                    'filename':         'xcf_exp_randomsonly.fits.gz',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_auto__lya_dla_cross/',
                    'result_name':      'result_lya_auto__lya_dla_cross_rmin40.0_rmax160.0_afree.h5',
                    'corr_name':        'LYA(LYA)xDLA',
                    #'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.6, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$\lambda_{r,\rm{DLA}}<1180 \AA$', 'leg_loc': 'shared'},
                    }
subplots[(1,1)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_056/picca_00296/',
                    'filename':         'xcf_exp_randomsonly.fits.gz',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_auto__lya_dla_cross/',
                    'result_name':      'result_lya_auto__lya_dla_cross_rmin40.0_rmax160.0_afree.h5',
                    'corr_name':        'LYA(LYA)xDLA',
                    #'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.6, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$\lambda_{r,\rm{DLA}}<1160 \AA$', 'leg_loc': 'shared'},
                    }
"""
"""
filename = 'lya_dla_cross_stack_contributions.pdf'
figsize=(18,8)
subplots = {}
subplots[(0,0)] =  {'location':         basedir+'/analysis/correlation_functions/stack/measurements/lya_dla_cross/',
                    'filename':         '/correlations/xcf_exp_lya_dla_cross.fits.gz',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_auto__lya_dla_cross/',
                    'result_name':      'result_lya_auto__lya_dla_cross_rmin40.0_rmax160.0_afree.h5',
                    'corr_name':        'LYA(LYA)xDLA',
                    #'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': 'standard removal', 'leg_loc': 'shared'},
                    }
subplots[(0,1)] =  {'location':         basedir+'/analysis/correlation_functions/stack/measurements/lya_dla_cross/',
                    'filename':         '/correlations/xcf_exp_lya_dla_cross_norandoms.fits.gz',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_auto__lya_dla_cross/',
                    'result_name':      'result_lya_auto__lya_dla_cross_rmin40.0_rmax160.0_afree.h5',
                    'corr_name':        'LYA(LYA)xDLA',
                    #'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': 'data only', 'leg_loc': 'shared'},
                    }
subplots[(0,2)] =  {'location':         basedir+'/analysis/correlation_functions/stack/measurements/lya_dla_cross/',
                    'filename':         '/correlations/xcf_exp_lya_dla_cross_randomsonly.fits.gz',
                    'result_location':  basedir+'/analysis/correlation_functions/stack/fits/lya_auto__lya_dla_cross/',
                    'result_name':      'result_lya_auto__lya_dla_cross_rmin40.0_rmax160.0_afree.h5',
                    'corr_name':        'LYA(LYA)xDLA',
                    #'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': 'randoms only', 'leg_loc': 'shared'},
                    }

"""
"""
#Create a dictionary with all information about the subplots:
filename = 'lya_dla_cross_compare_dlalrcuts.pdf'
figsize=(16,16)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_053/picca_00281/',
                    'filename':         'xcf_exp_randoms.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$\lambda_{r,DLA} < \lambda_\alpha$', 'leg_loc': 'shared'},
                    }
subplots[(0,1)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_054/picca_00286/',
                    'filename':         'xcf_exp_randoms.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$\lambda_{r,DLA} < 1200 \AA$', 'leg_loc': 'shared'},
                    }
subplots[(1,0)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_055/picca_00291/',
                    'filename':         'xcf_exp_randoms.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$\lambda_{r,DLA} < 1180 \AA$', 'leg_loc': 'shared'},
                    }
subplots[(1,1)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_056/picca_00296/',
                    'filename':         'xcf_exp_randoms.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$\lambda_{r,DLA} < 1160 \AA$', 'leg_loc': 'shared'},
                    }
"""
"""
filename = 'lya_dla_cross_compare_dlalrbins.pdf'
figsize=(16,16)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_057/picca_00301/',
                    'filename':         'xcf_exp_randoms.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$1180 < \lambda_{r,DLA}\ [\AA]\ < \lambda_\alpha$', 'leg_loc': 'shared'},
                    }
subplots[(0,1)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_058/picca_00306/',
                    'filename':         'xcf_exp_randoms.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$1140 < \lambda_{r,DLA}\ [\AA]\ < 1180$', 'leg_loc': 'shared'},
                    }
subplots[(1,0)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_059/picca_00311/',
                    'filename':         'xcf_exp_randoms.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$1100 < \lambda_{r,DLA}\ [\AA]\ < 1140$', 'leg_loc': 'shared'},
                    }
subplots[(1,1)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_060/picca_00316/',
                    'filename':         'xcf_exp_randoms.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$1060 < \lambda_{r,DLA}\ [\AA]\ < 1100$', 'leg_loc': 'shared'},
                    }
"""
"""
filename = 'lya_dla_cross_compare_deltalr1200_deltadlalr1200.pdf'
figsize=(16,8)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_050/picca_00266/',
                    'filename':         'xcf_exp_0.2_randoms.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$\lambda_r < 1200 \AA$', 'leg_loc': 'shared'},
                    }
subplots[(0,1)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_052/picca_00276/',
                    'filename':         'xcf_exp_0.2_randoms.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$\lambda_r < 1200 \AA,\ \lambda_{r,DLA} < 1200 \AA$', 'leg_loc': 'shared'},
                    }
"""
"""
filename = 'lya_dla_cross_zbins_compare_deltalr1200_deltadlalr1200.pdf'
figsize=(24,12)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_050/picca_00267/',
                    'filename':         'xcf_exp_0.2_randoms_0.0_2.2.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.09605, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$z_{\rm{QSO}}<2.2$', 'leg_loc': 'shared'},
                    }
subplots[(0,1)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_050/picca_00268/',
                    'filename':         'xcf_exp_0.2_randoms_2.2_2.6.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1161, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$2.2<z_{\rm{QSO}}<2.6$', 'leg_loc': 'shared'},
                    }
subplots[(0,2)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_050/picca_00269/',
                    'filename':         'xcf_exp_0.2_randoms_2.6_3.0.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1591, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$2.6<z_{\rm{QSO}}<3.0$', 'leg_loc': 'shared'},
                    }
subplots[(0,3)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_050/picca_00270/',
                    'filename':         'xcf_exp_0.2_randoms_3.0_10.0.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.2133, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$3.0<z_{\rm{QSO}}$', 'leg_loc': 'shared'},
                    }
subplots[(1,0)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_052/picca_00277/',
                    'filename':         'xcf_exp_0.2_randoms_0.0_2.2.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.09605, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$z_{\rm{QSO}}<2.2$', 'leg_loc': 'shared'},
                    }
subplots[(1,1)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_052/picca_00278/',
                    'filename':         'xcf_exp_0.2_randoms_2.2_2.6.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1161, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$2.2<z_{\rm{QSO}}<2.6$', 'leg_loc': 'shared'},
                    }
subplots[(1,2)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_052/picca_00279/',
                    'filename':         'xcf_exp_0.2_randoms_2.6_3.0.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1591, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$2.6<z_{\rm{QSO}}<3.0$', 'leg_loc': 'shared'},
                    }
subplots[(1,3)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_052/picca_00280/',
                    'filename':         'xcf_exp_0.2_randoms_3.0_10.0.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.2133, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$3.0<z_{\rm{QSO}}$', 'leg_loc': 'shared'},
                    }
"""
"""
filename = 'lya_dla_cross_zbins_new0.5.pdf'
figsize=(24,6)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_049/picca_00257/',
                    'filename':         'xcf_exp_randoms_0.0_2.2.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.09605, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$z_{\rm{QSO}}<2.2$', 'leg_loc': 'shared'},
                    }
subplots[(0,1)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_049/picca_00258/',
                    'filename':         'xcf_exp_randoms_2.2_2.6.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1161, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$2.2<z_{\rm{QSO}}<2.6$', 'leg_loc': 'shared'},
                    }
subplots[(0,2)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_049/picca_00259/',
                    'filename':         'xcf_exp_randoms_2.6_3.0.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1591, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$2.6<z_{\rm{QSO}}<3.0$', 'leg_loc': 'shared'},
                    }
subplots[(0,3)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_049/picca_00260/',
                    'filename':         'xcf_exp_randoms_3.0_10.0.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.2133, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$3.0<z_{\rm{QSO}}$', 'leg_loc': 'shared'},
                    }
"""
"""
filename = 'lya_dla_cross_zbins_new0.5_contributions.pdf'
figsize=(24,12)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_049/picca_00257/',
                    'filename':         'xcf_exp_norandoms_0.0_2.2.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.09605, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$z_{\rm{QSO}}<2.2$', 'leg_loc': 'shared'},
                    }
subplots[(0,1)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_049/picca_00258/',
                    'filename':         'xcf_exp_norandoms_2.2_2.6.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1161, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$2.2<z_{\rm{QSO}}<2.6$', 'leg_loc': 'shared'},
                    }
subplots[(0,2)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_049/picca_00259/',
                    'filename':         'xcf_exp_norandoms_2.6_3.0.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1591, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$2.6<z_{\rm{QSO}}<3.0$', 'leg_loc': 'shared'},
                    }
subplots[(0,3)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_049/picca_00260/',
                    'filename':         'xcf_exp_norandoms_3.0_10.0.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.2133, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$3.0<z_{\rm{QSO}}$', 'leg_loc': 'shared'},
                    }
subplots[(1,0)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_049/picca_00257/',
                    'filename':         'xcf_exp_randomsonly_0.0_2.2.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.09605, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$z_{\rm{QSO}}<2.2$', 'leg_loc': 'shared'},
                    }
subplots[(1,1)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_049/picca_00258/',
                    'filename':         'xcf_exp_randomsonly_2.2_2.6.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1161, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$2.2<z_{\rm{QSO}}<2.6$', 'leg_loc': 'shared'},
                    }
subplots[(1,2)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_049/picca_00259/',
                    'filename':         'xcf_exp_randomsonly_2.6_3.0.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1591, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$2.6<z_{\rm{QSO}}<3.0$', 'leg_loc': 'shared'},
                    }
subplots[(1,3)] =  {'location':         '/global/cscratch1/sd/jfarr/picca_output/picca_analysis_049/picca_00260/',
                    'filename':         'xcf_exp_randomsonly_3.0_10.0.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.2133, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': r'$3.0<z_{\rm{QSO}}$', 'leg_loc': 'shared'},
                    }
"""
"""
filename = 'dla_auto_vs_rp.pdf'
figsize=(12,8)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_046/picca_00251/',
                    'filename':         'co_exp.fits.gz',
                    'rt_bins':          [(0.0,16.0)],
                    #'rt_bins':          [(0.0,4.0),(16.0,20.0),(36.0,40.0),(56.0,60.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': 2.0, 'b2': 2.0, 'beta1': 0.48, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': '', 'leg_loc': 0},
                    }
"""
"""
filename = 'lya_dla_cross_vs_rp_dlalr.pdf'
figsize=(12,8)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_047/picca_00253/',
                    'filename':         'xcf_exp_0.2_randoms.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0),(12.0,16.0)],
                    #'rt_bins':          [(0.0,4.0),(16.0,20.0),(36.0,40.0),(56.0,60.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': '', 'leg_loc': 0},
                    }
"""
"""
filename = 'berkeley_cross_0.2_QSO_vs_rp.pdf'
figsize=(12,8)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_043/picca_00230/',
                    'filename':         'xcf_exp_0.2_shuffle.fits.gz',
                    #'rt_bins':          [(0.0,4.0),(48.0,52.0),(96.0,104.0),(148.0,152.0)],
                    'rt_bins':          [(0.0,4.0),(12.0,16.0),(24.0,28.0),(36.0,40.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 10., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': None},
                    }
"""
"""
filename = 'compare_xcf.pdf'
figsize=(12,6)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_044/picca_00238/',
                    'filename':         'xcf_exp_0.2_shuffle.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 10., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True},
                    }

subplots[(0,1)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_044/picca_00238_helion/',
                    'filename':         'xcf_z_0_10-exp.fits',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0)],
                    'ns':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   True,
                    'result_name':      'xcf_z_0_10-exp.h5',
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True},
                    }
"""
"""
filename = 'lya_dla_cross_contributions.pdf'
figsize=(12,6)
subplots = {}
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_045/picca_00240/',
                    'filename':         'xcf_exp_0.2_randoms.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': 'standard randoms removal', 'leg_loc': 0},
                    }
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_045/picca_00240/',
                    'filename':         'xcf_exp_0.2_norandoms.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': 'data only', 'leg_loc': 0},
                    }
subplots[(0,1)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_045/picca_00240/',
                    'filename':         'xcf_exp_0.2_randoms_only.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': 'randoms only', 'leg_loc': 0},
                    }
subplots[(0,1)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_045/picca_00240/',
                    'filename':         'xcf_exp_0.2_randoms_subbin.fits.gz',
                    'rt_bins':          [(0.0,4.0),(4.0,8.0),(8.0,12.0)],
                    #'rt_bins':          [(28.0,32.0),(32.0,36.0),(36.0,40.0),(40.0,44.0),(44.0,48.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3','C4'],
                    'plot_data':        {'r_power': 0, 'nr': 40},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'rmax': 160., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.119, 'b2': 2.0, 'beta1': 1.53, 'beta2': 0.48},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True, 'title': 'alt randoms removal', 'leg_loc': 0},
                    }
"""
"""
filename = 'corr_plot_vs_rt.pdf'
subplots = {}
subplots[(0,0)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_037/combined/',
                    'filename':         'cf_exp.fits.gz',
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0),(60.0,64.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax': 160.0},
                    'plot_picca_fit':   True,
                    'picca_fit_data':   {'rmin': 40., 'afix': 'free'},
                    'plot_manual_fit':  False,
                    'manual_fit_data':  {'b1': -0.133, 'b2': -0.133, 'beta1': 1.4, 'beta2': 1.4},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True},
                    }
subplots[(0,1)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_041/picca_00223/',
                    'filename':         'xcf_exp_2400k_shuffle.fits.gz',
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0),(60.0,64.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3'],
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
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0),(60.0,64.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax': 160.0},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1049, 'b2': -0.1049, 'beta1': 1.3783, 'beta2': 1.3783},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True},
                    }
subplots[(0,1)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_041/picca_00214/',
                    'filename':         'xcf_exp_noshuffle.fits.gz',
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0),(60.0,64.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3'],
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
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0),(60.0,64.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3'],
                    'plot_data':        {'r_power': 2, 'nr': 40, 'rmax': 160.0},
                    'plot_picca_fit':   False,
                    'picca_fit_data':   {'rmin': 40., 'afix': 'free'},
                    'plot_manual_fit':  True,
                    'manual_fit_data':  {'b1': -0.1049, 'b2': 2.0, 'beta1': 1.3783, 'beta2': 0.79},
                    'format':           {'legend': True, 'xlabel': True, 'ylabel': True},
                    }
subplots[(0,1)] =  {'location':         '/global/homes/j/jfarr/Programs/picca/picca_analysis_041/picca_00227/',
                    'filename':         'xcf_exp_zb0.05_2400k_pdcov.fits.gz',
                    'rt_bins':          [(0.0,4.0),(20.0,24.0),(40.0,44.0),(60.0,64.0)],
                    'rt_bin_colours':   ['C0','C1','C2','C3'],
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
    if subplots[key]['format']['leg_loc'] == 'shared':
        share_legend = True
    subplots[key]['corr_object'] = corr_obj
    plot_functions.plot_rt_bins_vs_rp(fig,axs[key],subplots[key],fontsize)
    #axs[key].set_ylim(-0.001,0.0025)
    axs[key].set_ylim(-0.004,0.004)
    axs[key].set_xlim(-200,200)
    #axs[key].grid()

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
