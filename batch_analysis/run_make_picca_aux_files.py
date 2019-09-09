#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import argparse
from scipy.interpolate import interp1d

from lyacolore import submit_utils

################################################################################

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#locations and names

parser.add_argument('--base-dir', type = str, default = None, required=True,
                    help = 'base directory')

parser.add_argument('--out-dir', type = str, default = None, required=False,
                    help = 'output directory')

parser.add_argument('--v-maj', type = int, default = 9, required=False,
                    help = 'major version of lyacolore realisations')

parser.add_argument('--v-min', type = int, default = 0, required=False,
                    help = 'minor version of lyacolore realisations')

parser.add_argument('--v-realisations', type = int, default = [0], required=False,
                    help = 'realisation numbers of lyacolore realisations', nargs='*')

parser.add_argument('--fid-Om', type=float, default=0.315, required=False,
                    help='Omega_matter(z=0) of fiducial LambdaCDM cosmology')

parser.add_argument('--fid-Or', type=float, default=0., required=False,
                    help='Omega_radiation(z=0) of fiducial LambdaCDM cosmology')

parser.add_argument('--fit-stack', action="store_true", default = False, required=False,
                    help = 'make files for fitting the stack of all realisations')

#Fit variables

parser.add_argument('--rmin-values', type = float, default = [20.], required=False,
                    help = 'minimum separation', nargs='*')

parser.add_argument('--rmax-values', type = float, default = [160.], required=False,
                    help = 'maximum separation', nargs='*')

parser.add_argument('--afix-values', type = str, default = ['free'], required=False,
                    help = 'maximum separation', nargs='*')

#Correlations to fit

parser.add_argument('--fit-lya-auto', action="store_true", default = False, required=False,
                    help = 'make files for fitting the lya auto correlation')

parser.add_argument('--fit-qso-auto', action="store_true", default = False, required=False,
                    help = 'make files for fitting the qso auto correlation')

parser.add_argument('--fit-dla-auto', action="store_true", default = False, required=False,
                    help = 'make files for fitting the dla auto correlation')

parser.add_argument('--fit-lya-aa-auto', action="store_true", default = False, required=False,
                    help = 'make files for fitting the lya all absorber auto correlation')

parser.add_argument('--fit-lya-qso-cross', action="store_true", default = False, required=False,
                    help = 'make files for fitting the lya qso cross correlation')

parser.add_argument('--fit-lya-dla-cross', action="store_true", default = False, required=False,
                    help = 'make files for fitting the lya dla cross correlation')

parser.add_argument('--fit-qso-dla-cross', action="store_true", default = False, required=False,
                    help = 'make files for fitting the qso dla cross correlation')

parser.add_argument('--fit-lya-auto--lya-qso-cross', action="store_true", default = False, required=False,
                    help = 'make files for joint lya-auto, lya-qso-cross fit')

parser.add_argument('--fit-lya-auto--lya-qso-cross--qso-auto', action="store_true", default = False, required=False,
                    help = 'make files for joint lya-auto, lya-qso-cross, qso-auto fit')

parser.add_argument('--fit-lya-auto--lya-dla-cross', action="store_true", default = False, required=False,
                    help = 'make files for joint lya-auto, lya-dla-cross fit')

parser.add_argument('--fit-lya-auto--lya-dla-cross--dla-auto', action="store_true", default = False, required=False,
                    help = 'make files for joint lya-auto, lya-dla-cross, dla-auto fit')

parser.add_argument('--fit-qso-auto--qso-dla-cross--dla-auto', action="store_true", default = False, required=False,
                    help = 'make files for joint qso-auto, qso-dla-cross, dla-auto fit')

parser.add_argument('--fit-all-correlations', action="store_true", default = False, required=False,
                    help = 'make files for the joint fit of all correlations')

parser.add_argument('--fit-all', action="store_true", default = False, required=False,
                    help = 'make files for fitting all correlations')

args = parser.parse_args()

################################################################################

if args.out_dir is None:
    args.out_dir = args.base_dir

if args.fit_all:
    args.fit_lya_auto = True
    args.fit_qso_auto = True
    args.fit_dla_auto = True
    args.fit_lya_aa_auto = True
    args.fit_lya_qso_cross = True
    args.fit_lya_dla_cross = True
    args.fit_qso_dla_cross = True
    args.fit_lya_auto__lya_qso_cross = True
    args.fit_lya_auto__lya_qso_cross__qso_auto = True
    args.fit_lya_auto__lya_dla_cross = True
    args.fit_lya_auto__lya_dla_cross__dla_auto = True
    args.fit_qso_auto__qso_dla_cross = True
    args.fit_qso_auto__qso_dla_cross__dla_auto = True
    args.fit_all_correlations = True

################################################################################

## Functions to make chi2, parameter and config files given the options.
def make_chi2_file(filepath,zeff,configs,result_filename):

    chi2_text = ''
    chi2_text += '#!/bin/bash -l\n\n'
    chi2_text += '[data sets]\n'
    chi2_text += 'zeff = {}\n'.format(zeff)
    config_string = ''
    for config in configs:
        config_string += config + ' '
    chi2_text += 'ini files = {}\n\n'.format(config_string)
    chi2_text += '[cosmo-fit type]\n'
    chi2_text += 'cosmo fit func = ap_at\n\n'
    chi2_text += '[verbosity]\n'
    chi2_text += 'level = 0\n\n'
    chi2_text += '[output]\n'
    chi2_text += 'filename = {}\n\n'.format(result_filename)
    chi2_text += '[fiducial]\n'
    chi2_text += 'filename = /PlanckDR12/PlanckDR12.fits\n'
    file = open(filepath,'w')
    file.write(chi2_text)
    file.close()

    return

def make_corr_info_file(filepath,corr_type,q1,q2):

    parameter_file_text = ''
    parameter_file_text += 'correl_type = {}\n'.format(corr_type)
    parameter_file_text += 'quantity 1 = {}\n'.format(q1)
    parameter_file_text += 'quantity 2 = {}\n'.format(q2)
    file = open(filepath,'w')
    file.write(parameter_file_text)
    file.close()

    return

def make_config_file(filepath,options_dict):

    config_text = ''
    config_text += '#!/bin/bash -l\n\n'

    for key in options_dict.keys():
        config_text += '[{}]\n'.format(key)
        sect_dict = options_dict[key]
        for option in sect_dict.keys():
            config_text += '{} = {}\n'.format(option,sect_dict[option])

    file = open(filepath,'w')
    file.write(config_text)
    file.close()

    return

## Functions to make chi2 and config files for each fit combination.
def make_lya_auto_fit_files(fits_dir,exp_filepaths,rmin=20.,rmax=160.,afix='free'):

    name = 'lya_auto'
    suffix = 'rmin{}_rmax{}_a{}'.format(rmin,rmax,afix)
    fit_dir = fits_dir + '/' + name + '/'
    submit_utils.check_dir(fit_dir)

    exp_files = [exp_filepaths[key] for key in exp_filepaths.keys() if key in name]
    zeff = get_zeff(exp_files)
    f = get_growth_rate(zeff,Om_z0=args.fid_Om)

    configs = []


    data_dict = {'name':            'LYA(LYA)xLYA(LYA)',
                 'tracer1':         'LYA',
                 'tracer2':         'LYA',
                 'tracer1-type':    'continuous',
                 'tracer2-type':    'continuous',
                 'filename':        exp_filepaths['lya_auto'],
                 'ell-max':         6,
                 }

    cuts_dict = {'rp-min':  0.,
                 'rp-max':  200.,
                 'rt-min':  0.,
                 'rt-max':  200.,
                 'r-min':   rmin,
                 'r-max':   rmax,
                 'mu-min':  0.,
                 'mu-max':  1.
                 }

    model_dict = {'model-pk':           'pk_kaiser',
                  'model-xi':           'xi',
                  'z evol LYA':         'bias_vs_z_std',
                  'growth function':    'growth_factor_de',
                  'pk-gauss-smoothing': 'pk_gauss_smoothing',
                  }

    parameters_dict = {'ap':                            '1.     0.1     0.8     1.2     {}'.format(afix),
                       'at':                            '1.     0.1     0.8     1.2     {}'.format(afix),
                       'bao_amp':                       '1.     0.      None    None    fixed',
                       'sigmaNL_per':                   '0.     0.      None    None    fixed',
                       'sigmaNL_par':                   '0.     0.      None    None    fixed',
                       'growth_rate':                   '{}     0.      None    None    fixed'.format(f),
                       'bias_eta_LYA':                  '-0.1   1.      None    None    free',
                       'beta_LYA':                      '0.5    1.      None    None    free',
                       'alpha_LYA':                     '2.9    0.      None    None    fixed',
                       'par binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                       'per binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                       'par_sigma_smooth':              '2.     2.      None    None    free',
                       'per_sigma_smooth':              '2.     2.      None    None    free',
                       }

    options_dict = {'data':         data_dict,
                    'cuts':         cuts_dict,
                    'model':        model_dict,
                    'parameters':   parameters_dict,
                    }

    config_filepath = fit_dir + '/config_lya_auto_{}.ini'.format(suffix)
    make_config_file(config_filepath,options_dict)
    configs += [config_filepath]

    chi2_filepath = fit_dir + 'chi2_{}_{}.ini'.format(name,suffix)
    result_filepath = fit_dir + 'result_{}_{}.h5'.format(name,suffix)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return options_dict

def make_qso_auto_fit_files(fits_dir,exp_filepaths,rmin=20.,rmax=160.,afix='free'):

    name = 'qso_auto'
    suffix = 'rmin{}_rmax{}_a{}'.format(rmin,rmax,afix)
    fit_dir = fits_dir + '/' + name + '/'
    submit_utils.check_dir(fit_dir)

    exp_files = [exp_filepaths[key] for key in exp_filepaths.keys() if key in name]
    zeff = get_zeff(exp_files)
    f = get_growth_rate(zeff,Om_z0=args.fid_Om)

    configs = []

    data_dict = {'name':            'QSOxQSO',
                 'tracer1':         'QSO',
                 'tracer2':         'QSO',
                 'tracer1-type':    'discrete',
                 'tracer2-type':    'discrete',
                 'filename':        exp_filepaths['qso_auto'],
                 'ell-max':         6,
                 }

    cuts_dict = {'rp-min':  0.,
                 'rp-max':  200.,
                 'rt-min':  0.,
                 'rt-max':  200.,
                 'r-min':   rmin,
                 'r-max':   rmax,
                 'mu-min':  0.,
                 'mu-max':  1.
                 }

    model_dict = {'model-pk':               'pk_kaiser',
                  'model-xi':               'xi_drp',
                  'z evol QSO':             'bias_vs_z_std',
                  'velocity dispersion':    'pk_velo_lorentz',
                  'growth function':        'growth_factor_de',
                  'pk-gauss-smoothing':     'pk_gauss_smoothing',
                  }

    parameters_dict = {'ap':                        '1.     0.1     0.8     1.2     {}'.format(afix),
                       'at':                        '1.     0.1     0.8     1.2     {}'.format(afix),
                       'bao_amp':                   '1.     0.      None    None    fixed',
                       'sigmaNL_per':               '0.     0.      None    None    fixed',
                       'sigmaNL_par':               '0.     0.      None    None    fixed',
                       'growth_rate':               '{}     0.      None    None    fixed'.format(f),
                       'bias_eta_QSO':              '1.     1.      None    None    fixed',
                       'beta_QSO':                  '0.25   1.      None    None    free',
                       'alpha_QSO':                 '1.44   0.      None    None    fixed',
                       'par binsize QSOxQSO':       '4      0.      None    None    fixed',
                       'per binsize QSOxQSO':       '4      0.      None    None    fixed',
                       'par_sigma_smooth':          '2.     2.      None    None    free',
                       'per_sigma_smooth':          '2.     2.      None    None    free',
                       'drp_QSO':                   '0.     0.1     None    None    fixed',
                       'sigma_velo_lorentz_QSO':    '0.     0.1     None    None    fixed',
                       }

    options_dict = {'data':         data_dict,
                    'cuts':         cuts_dict,
                    'model':        model_dict,
                    'parameters':   parameters_dict,
                    }

    config_filepath = fit_dir + '/config_qso_auto_{}.ini'.format(suffix)
    make_config_file(config_filepath,options_dict)
    configs += [config_filepath]

    chi2_filepath = fit_dir + 'chi2_{}_{}.ini'.format(name,suffix)
    result_filepath = fit_dir + 'result_{}_{}.h5'.format(name,suffix)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_dla_auto_fit_files(fits_dir,exp_filepaths,rmin=20.,rmax=160.,afix='free'):

    name = 'dla_auto'
    suffix = 'rmin{}_rmax{}_a{}'.format(rmin,rmax,afix)
    fit_dir = fits_dir + '/' + name + '/'
    submit_utils.check_dir(fit_dir)

    exp_files = [exp_filepaths[key] for key in exp_filepaths.keys() if key in name]
    zeff = get_zeff(exp_files)
    f = get_growth_rate(zeff,Om_z0=args.fid_Om)

    configs = []

    data_dict = {'name':            'DLAxDLA',
                 'tracer1':         'DLA',
                 'tracer2':         'DLA',
                 'tracer1-type':    'discrete',
                 'tracer2-type':    'discrete',
                 'filename':        exp_filepaths['dla_auto'],
                 'ell-max':         6,
                 }

    cuts_dict = {'rp-min':  0.,
                 'rp-max':  200.,
                 'rt-min':  0.,
                 'rt-max':  200.,
                 'r-min':   rmin,
                 'r-max':   rmax,
                 'mu-min':  0.,
                 'mu-max':  1.
                 }

    model_dict = {'model-pk':               'pk_kaiser',
                  'model-xi':               'xi_drp',
                  'z evol DLA':             'bias_vs_z_std',
                  'velocity dispersion':    'pk_velo_lorentz',
                  'growth function':        'growth_factor_de',
                  'pk-gauss-smoothing':     'pk_gauss_smoothing',
                  }

    parameters_dict = {'ap':                        '1.     0.1     0.8     1.2     {}'.format(afix),
                       'at':                        '1.     0.1     0.8     1.2     {}'.format(afix),
                       'bao_amp':                   '1.     0.      None    None    fixed',
                       'sigmaNL_per':               '0.     0.      None    None    fixed',
                       'sigmaNL_par':               '0.     0.      None    None    fixed',
                       'growth_rate':               '{}     0.      None    None    fixed'.format(f),
                       'bias_eta_DLA':              '1.     1.      None    None    fixed',
                       'beta_DLA':                  '0.5    1.      None    None    free',
                       'alpha_DLA':                 '0.0    0.      None    None    fixed',
                       'par binsize DLAxDLA':       '4      0.      None    None    fixed',
                       'per binsize DLAxDLA':       '4      0.      None    None    fixed',
                       'par_sigma_smooth':          '2.     2.      None    None    free',
                       'per_sigma_smooth':          '2.     2.      None    None    free',
                       'drp_DLA':                   '0.     0.1     None    None    fixed',
                       'sigma_velo_lorentz_DLA':    '0.     0.1     None    None    fixed',
                       }

    options_dict = {'data':         data_dict,
                    'cuts':         cuts_dict,
                    'model':        model_dict,
                    'parameters':   parameters_dict,
                    }

    config_filepath = fit_dir + '/config_dla_auto_{}.ini'.format(suffix)
    make_config_file(config_filepath,options_dict)
    configs += [config_filepath]

    chi2_filepath = fit_dir + 'chi2_{}_{}.ini'.format(name,suffix)
    result_filepath = fit_dir + 'result_{}_{}.h5'.format(name,suffix)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_lya_aa_auto_fit_files(fits_dir,exp_filepaths,rmin=20.,rmax=160.,afix='free'):

    name = 'lya_aa_auto'
    suffix = 'rmin{}_rmax{}_a{}'.format(rmin,rmax,afix)
    fit_dir = fits_dir + '/' + name + '/'
    submit_utils.check_dir(fit_dir)

    exp_files = [exp_filepaths[key] for key in exp_filepaths.keys() if key in name]
    zeff = get_zeff(exp_files)
    f = get_growth_rate(zeff,Om_z0=args.fid_Om)

    configs = []

    data_dict = {'name':            'LYA(LYA)xLYA(LYA)',
                 'tracer1':         'LYA',
                 'tracer2':         'LYA',
                 'tracer1-type':    'continuous',
                 'tracer2-type':    'continuous',
                 'filename':        exp_filepaths['lya_aa_auto'],
                 'ell-max':         6,
                 }

    cuts_dict = {'rp-min':  0.,
                 'rp-max':  200.,
                 'rt-min':  0.,
                 'rt-max':  200.,
                 'r-min':   rmin,
                 'r-max':   rmax,
                 'mu-min':  0.,
                 'mu-max':  1.
                 }

    model_dict = {'model-pk':           'pk_kaiser',
                  'model-xi':           'xi',
                  'z evol LYA':         'bias_vs_z_std',
                  'growth function':    'growth_factor_de',
                  'pk-gauss-smoothing': 'pk_gauss_smoothing',
                  }

    parameters_dict = {'ap':                            '1.     0.1     0.8     1.2     {}'.format(afix),
                       'at':                            '1.     0.1     0.8     1.2     {}'.format(afix),
                       'bao_amp':                       '1.     0.      None    None    fixed',
                       'sigmaNL_per':                   '0.     0.      None    None    fixed',
                       'sigmaNL_par':                   '0.     0.      None    None    fixed',
                       'growth_rate':                   '{}     0.      None    None    fixed'.format(f),
                       'bias_eta_LYA':                  '-0.1   1.      None    None    free',
                       'beta_LYA':                      '0.5    1.      None    None    free',
                       'alpha_LYA':                     '2.9    0.      None    None    fixed',
                       'par binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                       'per binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                       'par_sigma_smooth':              '2.     2.      None    None    free',
                       'per_sigma_smooth':              '2.     2.      None    None    free',
                       'bias_eta_SiII(1260)':           '-0.47E-3   0.01    None    None    free',
                       'beta_SiII(1260)':               '0.5        0.      None    None    fixed',
                       'alpha_SiII(1260)':              '1.0        0.      None    None    fixed',
                       'bias_eta_SiIII(1207)':          '-2.02E-3   0.01    None    None    free',
                       'beta_SiIII(1207)':              '0.5        0.      None    None    fixed',
                       'alpha_SiIII(1207)':             '1.0        0.      None    None    fixed',
                       'bias_eta_SiII(1193)':           '-1.18E-3   0.01    None    None    free',
                       'beta_SiII(1193)':               '0.5        0.      None    None    fixed',
                       'alpha_SiII(1193)':              '1.0        0.      None    None    fixed',
                       'bias_eta_SiII(1190)':           '-0.47E-3   0.01    None    None    free',
                       'beta_SiII(1190)':               '0.5        0.      None    None    fixed',
                       'alpha_SiII(1190)':              '1.0        0.      None    None    fixed',
                       'bias_eta_LYB':                  '-0.47E-3   0.01    None    None    free',
                       'beta_LYB':                      '0.5        0.      None    None    fixed',
                       'alpha_LYB':                     '1.0        0.      None    None    fixed',
                       }

    metals_dict = {'filename':      '{}/data/additional_data/metal_dmat_z_0_10.fits.gz\n'.format(args.base_dir),
                   'model-pk-met':  'pk_kaiser',
                   'model-xi-met':  'cached_xi_kaiser',
                   'z evol':        'bias_vs_z_std',
                   'in tracer1':    'SiII(1260) SiIII(1207) SiII(1193) SiII(1190)',
                   'in tracer2':    'SiII(1260) SiIII(1207) SiII(1193) SiII(1190)',
                   }

    options_dict = {'data':         data_dict,
                    'cuts':         cuts_dict,
                    'model':        model_dict,
                    'parameters':   parameters_dict,
                    'metals':       metals_dict,
                    }

    config_filepath = fit_dir + '/config_lya_aa_auto_{}.ini'.format(suffix)
    make_config_file(config_filepath,options_dict)
    configs += [config_filepath]

    chi2_filepath = fit_dir + 'chi2_{}_{}.ini'.format(name,suffix)
    result_filepath = fit_dir + 'result_{}_{}.h5'.format(name,suffix)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_lya_qso_cross_fit_files(fits_dir,exp_filepaths,rmin=20.,rmax=160.,afix='free'):

    name = 'lya_qso_cross'
    suffix = 'rmin{}_rmax{}_a{}'.format(rmin,rmax,afix)
    fit_dir = fits_dir + '/' + name + '/'
    submit_utils.check_dir(fit_dir)

    exp_files = [exp_filepaths[key] for key in exp_filepaths.keys() if key in name]
    zeff = get_zeff(exp_files)
    f = get_growth_rate(zeff,Om_z0=args.fid_Om)
    beta_QSO = get_beta_obj(zeff,obj='QSO')

    configs = []

    data_dict = {'name':            'LYA(LYA)xQSO',
                 'tracer1':         'LYA',
                 'tracer2':         'QSO',
                 'tracer1-type':    'continuous',
                 'tracer2-type':    'discrete',
                 'filename':        exp_filepaths['lya_qso_cross'],
                 'ell-max':         6,
                 }

    cuts_dict = {'rp-min':  -200.,
                 'rp-max':  200.,
                 'rt-min':  0.,
                 'rt-max':  200.,
                 'r-min':   rmin,
                 'r-max':   rmax,
                 'mu-min':  -1.,
                 'mu-max':  1.
                 }

    model_dict = {'model-pk':               'pk_kaiser',
                  'model-xi':               'xi_drp',
                  'z evol LYA':             'bias_vs_z_std',
                  'velocity dispersion':    'pk_velo_lorentz',
                  'z evol QSO':             'bias_vs_z_std',
                  'growth function':        'growth_factor_de',
                  'pk-gauss-smoothing':     'pk_gauss_smoothing',
                  }

    parameters_dict = {'ap':                        '1.     0.1     0.8     1.2     {}'.format(afix),
                       'at':                        '1.     0.1     0.8     1.2     {}'.format(afix),
                       'bao_amp':                   '1.     0.      None    None    fixed',
                       'sigmaNL_per':               '0.     0.      None    None    fixed',
                       'sigmaNL_par':               '0.     0.      None    None    fixed',
                       'growth_rate':               '{}     0.      None    None    fixed'.format(f),
                       'bias_eta_LYA':              '-0.1   1.      None    None    free',
                       'beta_LYA':                  '0.5    1.      None    None    free',
                       'alpha_LYA':                 '2.9    0.      None    None    fixed',
                       'bias_eta_QSO':              '1.     1.      None    None    fixed',
                       'beta_QSO':                  '{}     1.      None    None    fixed'.format(beta_QSO), # TODO: NEED TO calculate
                       'alpha_QSO':                 '1.44   0.      None    None    fixed',
                       'par binsize LYA(LYA)xQSO':  '4      0.      None    None    fixed',
                       'per binsize LYA(LYA)xQSO':  '4      0.      None    None    fixed',
                       'par_sigma_smooth':          '2.     2.      None    None    free',
                       'per_sigma_smooth':          '2.     2.      None    None    free',
                       'drp_QSO':                   '0.     0.1     None    None    fixed',
                       'sigma_velo_lorentz_QSO':    '0.     0.1     None    None    fixed',
                       }

    options_dict = {'data':         data_dict,
                    'cuts':         cuts_dict,
                    'model':        model_dict,
                    'parameters':   parameters_dict,
                    }

    config_filepath = fit_dir + '/config_lya_qso_cross_{}.ini'.format(suffix)
    make_config_file(config_filepath,options_dict)
    configs += [config_filepath]

    chi2_filepath = fit_dir + 'chi2_{}_{}.ini'.format(name,suffix)
    result_filepath = fit_dir + 'result_{}_{}.h5'.format(name,suffix)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_lya_dla_cross_fit_files(fits_dir,exp_filepaths,rmin=20.,rmax=160.,afix='free'):

    name = 'lya_dla_cross'
    suffix = 'rmin{}_rmax{}_a{}'.format(rmin,rmax,afix)
    fit_dir = fits_dir + '/' + name + '/'
    submit_utils.check_dir(fit_dir)

    exp_files = [exp_filepaths[key] for key in exp_filepaths.keys() if key in name]
    zeff = get_zeff(exp_files)
    f = get_growth_rate(zeff,Om_z0=args.fid_Om)
    beta_DLA = get_beta_obj(zeff,obj='DLA')

    configs = []

    data_dict = {'name':            'LYA(LYA)xDLA',
                 'tracer1':         'LYA',
                 'tracer2':         'DLA',
                 'tracer1-type':    'continuous',
                 'tracer2-type':    'discrete',
                 'filename':        exp_filepaths['lya_dla_cross'],
                 'ell-max':         6,
                 }

    cuts_dict = {'rp-min':  -200.,
                 'rp-max':  200.,
                 'rt-min':  0.,
                 'rt-max':  200.,
                 'r-min':   rmin,
                 'r-max':   rmax,
                 'mu-min':  -1.,
                 'mu-max':  1.
                 }

    model_dict = {'model-pk':               'pk_kaiser',
                  'model-xi':               'xi_drp',
                  'z evol LYA':             'bias_vs_z_std',
                  'velocity dispersion':    'pk_velo_lorentz',
                  'z evol DLA':             'bias_vs_z_std',
                  'growth function':        'growth_factor_de',
                  'pk-gauss-smoothing':     'pk_gauss_smoothing',
                  }

    parameters_dict = {'ap':                        '1.     0.1     0.8     1.2     {}'.format(afix),
                       'at':                        '1.     0.1     0.8     1.2     {}'.format(afix),
                       'bao_amp':                   '1.     0.      None    None    fixed',
                       'sigmaNL_per':               '0.     0.      None    None    fixed',
                       'sigmaNL_par':               '0.     0.      None    None    fixed',
                       'growth_rate':               '{}     0.      None    None    fixed'.format(f),
                       'bias_eta_LYA':              '-0.1   1.      None    None    free',
                       'beta_LYA':                  '0.5    1.      None    None    free',
                       'alpha_LYA':                 '2.9    0.      None    None    fixed',
                       'bias_eta_DLA':              '1.     1.      None    None    fixed',
                       'beta_DLA':                  '{}     1.      None    None    fixed'.format(beta_DLA),
                       'alpha_DLA':                 '0.0    0.      None    None    fixed',
                       'par binsize LYA(LYA)xDLA':  '4      0.      None    None    fixed',
                       'per binsize LYA(LYA)xDLA':  '4      0.      None    None    fixed',
                       'par_sigma_smooth':          '2.     2.      None    None    free',
                       'per_sigma_smooth':          '2.     2.      None    None    free',
                       'drp_DLA':                   '0.     0.1     None    None    fixed',
                       'sigma_velo_lorentz_DLA':    '0.     0.1     None    None    fixed',
                       }

    options_dict = {'data':         data_dict,
                    'cuts':         cuts_dict,
                    'model':        model_dict,
                    'parameters':   parameters_dict,
                    }

    config_filepath = fit_dir + '/config_lya_dla_cross_{}.ini'.format(suffix)
    make_config_file(config_filepath,options_dict)
    configs += [config_filepath]

    chi2_filepath = fit_dir + 'chi2_{}_{}.ini'.format(name,suffix)
    result_filepath = fit_dir + 'result_{}_{}.h5'.format(name,suffix)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_qso_dla_cross_fit_files(fits_dir,exp_filepaths,rmin=20.,rmax=160.,afix='free'):

    name = 'qso_dla_cross'
    suffix = 'rmin{}_rmax{}_a{}'.format(rmin,rmax,afix)
    fit_dir = fits_dir + '/' + name + '/'
    submit_utils.check_dir(fit_dir)

    exp_files = [exp_filepaths[key] for key in exp_filepaths.keys() if key in name]
    zeff = get_zeff(exp_files)
    f = get_growth_rate(zeff,Om_z0=args.fid_Om)

    configs = []

    data_dict = {'name':            'QSOxDLA',
                 'tracer1':         'QSO',
                 'tracer2':         'DLA',
                 'tracer1-type':    'continuous',
                 'tracer2-type':    'discrete',
                 'filename':        exp_filepaths['qso_dla_cross'],
                 'ell-max':         6,
                 }

    cuts_dict = {'rp-min':  -200.,
                 'rp-max':  200.,
                 'rt-min':  0.,
                 'rt-max':  200.,
                 'r-min':   rmin,
                 'r-max':   rmax,
                 'mu-min':  -1.,
                 'mu-max':  1.
                 }

    model_dict = {'model-pk':               'pk_kaiser',
                  'model-xi':               'xi_drp',
                  'z evol QSO':             'bias_vs_z_std',
                  'velocity dispersion':    'pk_velo_lorentz',
                  'z evol DLA':             'bias_vs_z_std',
                  'growth function':        'growth_factor_de',
                  'pk-gauss-smoothing':     'pk_gauss_smoothing',
                  }

    parameters_dict = {'ap':                        '1.     0.1     0.8     1.2     {}'.format(afix),
                       'at':                        '1.     0.1     0.8     1.2     {}'.format(afix),
                       'bao_amp':                   '1.     0.      None    None    fixed',
                       'sigmaNL_per':               '0.     0.      None    None    fixed',
                       'sigmaNL_par':               '0.     0.      None    None    fixed',
                       'growth_rate':               '{}     0.      None    None    fixed'.format(f),
                       'bias_eta_QSO':              '1.     1.      None    None    fixed',
                       'beta_QSO':                  '0.5    1.      None    None    free', # TODO: can we leave both betas free?
                       'alpha_QSO':                 '1.44   0.      None    None    fixed',
                       'bias_eta_DLA':              '1.     1.      None    None    fixed',
                       'beta_DLA':                  '0.25   1.      None    None    free', # TODO: can we leave both betas free?
                       'alpha_DLA':                 '0.0    0.      None    None    fixed',
                       'par binsize QSOxDLA':       '4      0.      None    None    fixed',
                       'per binsize QSOxDLA':       '4      0.      None    None    fixed',
                       'par_sigma_smooth':          '2.     2.      None    None    free',
                       'per_sigma_smooth':          '2.     2.      None    None    free',
                       'drp_QSO':                   '0.     0.1     None    None    fixed', # TODO: can we leave both drps free?
                       'sigma_velo_lorentz_QSO':    '0.     0.1     None    None    fixed',
                       'drp_DLA':                   '0.     0.1     None    None    fixed', # TODO: can we leave both drps free?
                       'sigma_velo_lorentz_DLA':    '0.     0.1     None    None    fixed',
                       }

    options_dict = {'data':         data_dict,
                    'cuts':         cuts_dict,
                    'model':        model_dict,
                    'parameters':   parameters_dict,
                    }

    config_filepath = fit_dir + '/config_qso_dla_cross_{}.ini'.format(suffix)
    make_config_file(config_filepath,options_dict)
    configs += [config_filepath]

    chi2_filepath = fit_dir + 'chi2_{}_{}.ini'.format(name,suffix)
    result_filepath = fit_dir + 'result_{}_{}.h5'.format(name,suffix)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_lya_auto__lya_qso_cross_fit_files(fits_dir,exp_filepaths,rmin=20.,rmax=160.,afix='free'):

    name = 'lya_auto__lya_qso_cross'
    suffix = 'rmin{}_rmax{}_a{}'.format(rmin,rmax,afix)
    fit_dir = fits_dir + '/' + name + '/'
    submit_utils.check_dir(fit_dir)

    exp_files = [exp_filepaths[key] for key in exp_filepaths.keys() if key in name]
    zeff = get_zeff(exp_files)
    f = get_growth_rate(zeff,Om_z0=args.fid_Om)

    configs = []

    lya_auto_data_dict = {'name':            'LYA(LYA)xLYA(LYA)',
                          'tracer1':         'LYA',
                          'tracer2':         'LYA',
                          'tracer1-type':    'continuous',
                          'tracer2-type':    'continuous',
                          'filename':        exp_filepaths['lya_auto'],
                          'ell-max':         6,
                          }

    lya_auto_cuts_dict = {'rp-min':  0.,
                          'rp-max':  200.,
                          'rt-min':  0.,
                          'rt-max':  200.,
                          'r-min':   rmin,
                          'r-max':   rmax,
                          'mu-min':  0.,
                          'mu-max':  1.
                          }

    lya_auto_model_dict = {'model-pk':           'pk_kaiser',
                           'model-xi':           'xi',
                           'z evol LYA':         'bias_vs_z_std',
                           'growth function':    'growth_factor_de',
                           'pk-gauss-smoothing': 'pk_gauss_smoothing',
                           }

    lya_auto_parameters_dict = {'ap':                            '1.     0.1     0.8     1.2     {}'.format(afix),
                                'at':                            '1.     0.1     0.8     1.2     {}'.format(afix),
                                'bao_amp':                       '1.     0.      None    None    fixed',
                                'sigmaNL_per':                   '0.     0.      None    None    fixed',
                                'sigmaNL_par':                   '0.     0.      None    None    fixed',
                                'growth_rate':                   '{}     0.      None    None    fixed'.format(f),
                                'bias_eta_LYA':                  '-0.1   1.      None    None    free',
                                'beta_LYA':                      '0.5    1.      None    None    free',
                                'alpha_LYA':                     '2.9    0.      None    None    fixed',
                                'par binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'per binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'par_sigma_smooth':              '2.     2.      None    None    free',
                                'per_sigma_smooth':              '2.     2.      None    None    free',
                                }

    lya_auto_options_dict = {'data':         lya_auto_data_dict,
                             'cuts':         lya_auto_cuts_dict,
                             'model':        lya_auto_model_dict,
                             'parameters':   lya_auto_parameters_dict,
                             }

    config_filepath = fit_dir + '/config_lya_auto_{}.ini'.format(suffix)
    make_config_file(config_filepath,lya_auto_options_dict)
    configs += [config_filepath]

    lya_qso_cross_data_dict = {'name':            'LYA(LYA)xQSO',
                               'tracer1':         'LYA',
                               'tracer2':         'QSO',
                               'tracer1-type':    'continuous',
                               'tracer2-type':    'discrete',
                               'filename':        exp_filepaths['lya_qso_cross'],
                               'ell-max':         6,
                               }

    lya_qso_cross_cuts_dict = {'rp-min':  -200.,
                               'rp-max':  200.,
                               'rt-min':  0.,
                               'rt-max':  200.,
                               'r-min':   rmin,
                               'r-max':   rmax,
                               'mu-min':  -1.,
                               'mu-max':  1.
                               }

    lya_qso_cross_model_dict = {'model-pk':               'pk_kaiser',
                                'model-xi':               'xi_drp',
                                'z evol LYA':             'bias_vs_z_std',
                                'velocity dispersion':    'pk_velo_lorentz',
                                'z evol QSO':             'bias_vs_z_std',
                                'growth function':        'growth_factor_de',
                                'pk-gauss-smoothing':     'pk_gauss_smoothing',
                                }

    lya_qso_cross_parameters_dict = {'bias_eta_QSO':              '1.     1.      None    None    fixed',
                                     'beta_QSO':                  '0.25   1.      None    None    free',
                                     'alpha_QSO':                 '1.44   0.      None    None    fixed',
                                     'par binsize LYA(LYA)xQSO':  '4      0.      None    None    fixed',
                                     'per binsize LYA(LYA)xQSO':  '4      0.      None    None    fixed',
                                     'drp_QSO':                   '0.     0.1     None    None    fixed',
                                     'sigma_velo_lorentz_QSO':    '0.     0.1     None    None    fixed',
                                     }

    lya_qso_cross_options_dict = {'data':         lya_qso_cross_data_dict,
                                  'cuts':         lya_qso_cross_cuts_dict,
                                  'model':        lya_qso_cross_model_dict,
                                  'parameters':   lya_qso_cross_parameters_dict,
                                  }

    config_filepath = fit_dir + '/config_lya_qso_cross_{}.ini'.format(suffix)
    make_config_file(config_filepath,lya_qso_cross_options_dict)
    configs += [config_filepath]

    chi2_filepath = fit_dir + 'chi2_{}_{}.ini'.format(name,suffix)
    result_filepath = fit_dir + 'result_{}_{}.h5'.format(name,suffix)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_lya_auto__lya_qso_cross__qso_auto_fit_files(fits_dir,exp_filepaths,rmin=20.,rmax=160.,afix='free'):

    name = 'lya_auto__lya_qso_cross__qso_auto'
    suffix = 'rmin{}_rmax{}_a{}'.format(rmin,rmax,afix)
    fit_dir = fits_dir + '/' + name + '/'
    submit_utils.check_dir(fit_dir)

    exp_files = [exp_filepaths[key] for key in exp_filepaths.keys() if key in name]
    zeff = get_zeff(exp_files)
    f = get_growth_rate(zeff,Om_z0=args.fid_Om)

    configs = []

    lya_auto_data_dict = {'name':            'LYA(LYA)xLYA(LYA)',
                          'tracer1':         'LYA',
                          'tracer2':         'LYA',
                          'tracer1-type':    'continuous',
                          'tracer2-type':    'continuous',
                          'filename':        exp_filepaths['lya_auto'],
                          'ell-max':         6,
                          }

    lya_auto_cuts_dict = {'rp-min':  0.,
                          'rp-max':  200.,
                          'rt-min':  0.,
                          'rt-max':  200.,
                          'r-min':   rmin,
                          'r-max':   rmax,
                          'mu-min':  0.,
                          'mu-max':  1.
                          }

    lya_auto_model_dict = {'model-pk':           'pk_kaiser',
                           'model-xi':           'xi',
                           'z evol LYA':         'bias_vs_z_std',
                           'growth function':    'growth_factor_de',
                           'pk-gauss-smoothing': 'pk_gauss_smoothing',
                           }

    lya_auto_parameters_dict = {'ap':                            '1.     0.1     0.8     1.2     {}'.format(afix),
                                'at':                            '1.     0.1     0.8     1.2     {}'.format(afix),
                                'bao_amp':                       '1.     0.      None    None    fixed',
                                'sigmaNL_per':                   '0.     0.      None    None    fixed',
                                'sigmaNL_par':                   '0.     0.      None    None    fixed',
                                'growth_rate':                   '{}     0.      None    None    free'.format(f),
                                'bias_eta_LYA':                  '-0.1   1.      None    None    free',
                                'beta_LYA':                      '0.5    1.      None    None    free',
                                'alpha_LYA':                     '2.9    0.      None    None    fixed',
                                'par binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'per binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'par_sigma_smooth':              '2.     2.      None    None    free',
                                'per_sigma_smooth':              '2.     2.      None    None    free',
                                }

    lya_auto_options_dict = {'data':         lya_auto_data_dict,
                             'cuts':         lya_auto_cuts_dict,
                             'model':        lya_auto_model_dict,
                             'parameters':   lya_auto_parameters_dict,
                             }

    config_filepath = fit_dir + '/config_lya_auto_{}.ini'.format(suffix)
    make_config_file(config_filepath,lya_auto_options_dict)
    configs += [config_filepath]

    lya_qso_cross_data_dict = {'name':            'LYA(LYA)xQSO',
                               'tracer1':         'LYA',
                               'tracer2':         'QSO',
                               'tracer1-type':    'continuous',
                               'tracer2-type':    'discrete',
                               'filename':        exp_filepaths['lya_qso_cross'],
                               'ell-max':         6,
                               }

    lya_qso_cross_cuts_dict = {'rp-min':  -200.,
                               'rp-max':  200.,
                               'rt-min':  0.,
                               'rt-max':  200.,
                               'r-min':   rmin,
                               'r-max':   rmax,
                               'mu-min':  -1.,
                               'mu-max':  1.
                               }

    lya_qso_cross_model_dict = {'model-pk':               'pk_kaiser',
                                'model-xi':               'xi_drp',
                                'z evol LYA':             'bias_vs_z_std',
                                'velocity dispersion':    'pk_velo_lorentz',
                                'z evol QSO':             'bias_vs_z_std',
                                'growth function':        'growth_factor_de',
                                'pk-gauss-smoothing':     'pk_gauss_smoothing',
                                }

    lya_qso_cross_parameters_dict = {'bias_eta_QSO':              '1.     1.      None    None    fixed',
                                     'beta_QSO':                  '0.25   1.      None    None    free',
                                     'alpha_QSO':                 '1.44   0.      None    None    fixed',
                                     'par binsize LYA(LYA)xQSO':  '4      0.      None    None    fixed',
                                     'per binsize LYA(LYA)xQSO':  '4      0.      None    None    fixed',
                                     'drp_QSO':                   '0.     0.1     None    None    fixed',
                                     'sigma_velo_lorentz_QSO':    '0.     0.1     None    None    fixed',
                                     }

    lya_qso_cross_options_dict = {'data':         lya_qso_cross_data_dict,
                                  'cuts':         lya_qso_cross_cuts_dict,
                                  'model':        lya_qso_cross_model_dict,
                                  'parameters':   lya_qso_cross_parameters_dict,
                                  }

    config_filepath = fit_dir + '/config_lya_qso_cross_{}.ini'.format(suffix)
    make_config_file(config_filepath,lya_qso_cross_options_dict)
    configs += [config_filepath]

    qso_auto_data_dict = {'name':            'QSOxQSO',
                          'tracer1':         'QSO',
                          'tracer2':         'QSO',
                          'tracer1-type':    'discrete',
                          'tracer2-type':    'discrete',
                          'filename':        exp_filepaths['qso_auto'],
                          'ell-max':         6,
                          }

    qso_auto_cuts_dict = {'rp-min':  0.,
                          'rp-max':  200.,
                          'rt-min':  0.,
                          'rt-max':  200.,
                          'r-min':   rmin,
                          'r-max':   rmax,
                          'mu-min':  0.,
                          'mu-max':  1.
                          }

    qso_auto_model_dict = {'model-pk':               'pk_kaiser',
                           'model-xi':               'xi_drp',
                           'velocity dispersion':    'pk_velo_lorentz',
                           'z evol QSO':             'bias_vs_z_std',
                           'growth function':        'growth_factor_de',
                           'pk-gauss-smoothing':     'pk_gauss_smoothing',
                           }

    qso_auto_parameters_dict = {'par binsize QSOxQSO':  '4      0.      None    None    fixed',
                                'per binsize QSOxQSO':  '4      0.      None    None    fixed',
                                }

    qso_auto_options_dict = {'data':         qso_auto_data_dict,
                             'cuts':         qso_auto_cuts_dict,
                             'model':        qso_auto_model_dict,
                             'parameters':   qso_auto_parameters_dict,
                             }

    config_filepath = fit_dir + '/config_qso_auto_{}.ini'.format(suffix)
    make_config_file(config_filepath,qso_auto_options_dict)
    configs += [config_filepath]

    chi2_filepath = fit_dir + 'chi2_{}_{}.ini'.format(name,suffix)
    result_filepath = fit_dir + 'result_{}_{}.h5'.format(name,suffix)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_lya_auto__lya_dla_cross_fit_files(fits_dir,exp_filepaths,rmin=20.,rmax=160.,afix='free'):

    name = 'lya_auto__lya_dla_cross'
    suffix = 'rmin{}_rmax{}_a{}'.format(rmin,rmax,afix)
    fit_dir = fits_dir + '/' + name + '/'
    submit_utils.check_dir(fit_dir)

    exp_files = [exp_filepaths[key] for key in exp_filepaths.keys() if key in name]
    zeff = get_zeff(exp_files)
    f = get_growth_rate(zeff,Om_z0=args.fid_Om)

    configs = []

    lya_auto_data_dict = {'name':            'LYA(LYA)xLYA(LYA)',
                          'tracer1':         'LYA',
                          'tracer2':         'LYA',
                          'tracer1-type':    'continuous',
                          'tracer2-type':    'continuous',
                          'filename':        exp_filepaths['lya_auto'],
                          'ell-max':         6,
                          }

    lya_auto_cuts_dict = {'rp-min':  0.,
                          'rp-max':  200.,
                          'rt-min':  0.,
                          'rt-max':  200.,
                          'r-min':   rmin,
                          'r-max':   rmax,
                          'mu-min':  0.,
                          'mu-max':  1.
                          }

    lya_auto_model_dict = {'model-pk':           'pk_kaiser',
                           'model-xi':           'xi',
                           'z evol LYA':         'bias_vs_z_std',
                           'growth function':    'growth_factor_de',
                           'pk-gauss-smoothing': 'pk_gauss_smoothing',
                           }

    lya_auto_parameters_dict = {'ap':                            '1.     0.1     0.8     1.2     {}'.format(afix),
                                'at':                            '1.     0.1     0.8     1.2     {}'.format(afix),
                                'bao_amp':                       '1.     0.      None    None    fixed',
                                'sigmaNL_per':                   '0.     0.      None    None    fixed',
                                'sigmaNL_par':                   '0.     0.      None    None    fixed',
                                'growth_rate':                   '{}     0.      None    None    fixed'.format(f),
                                'bias_eta_LYA':                  '-0.1   1.      None    None    free',
                                'beta_LYA':                      '0.5    1.      None    None    free',
                                'alpha_LYA':                     '2.9    0.      None    None    fixed',
                                'par binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'per binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'par_sigma_smooth':              '2.     2.      None    None    free',
                                'per_sigma_smooth':              '2.     2.      None    None    free',
                                }

    lya_auto_options_dict = {'data':         lya_auto_data_dict,
                             'cuts':         lya_auto_cuts_dict,
                             'model':        lya_auto_model_dict,
                             'parameters':   lya_auto_parameters_dict,
                             }

    config_filepath = fit_dir + '/config_lya_auto_{}.ini'.format(suffix)
    make_config_file(config_filepath,lya_auto_options_dict)
    configs += [config_filepath]

    lya_dla_cross_data_dict = {'name':            'LYA(LYA)xDLA',
                               'tracer1':         'LYA',
                               'tracer2':         'DLA',
                               'tracer1-type':    'continuous',
                               'tracer2-type':    'discrete',
                               'filename':        exp_filepaths['lya_dla_cross'],
                               'ell-max':         6,
                               }

    lya_dla_cross_cuts_dict = {'rp-min':  -200.,
                               'rp-max':  200.,
                               'rt-min':  0.,
                               'rt-max':  200.,
                               'r-min':   rmin,
                               'r-max':   rmax,
                               'mu-min':  -1.,
                               'mu-max':  1.
                               }

    lya_dla_cross_model_dict = {'model-pk':               'pk_kaiser',
                                'model-xi':               'xi_drp',
                                'z evol LYA':             'bias_vs_z_std',
                                'velocity dispersion':    'pk_velo_lorentz',
                                'z evol DLA':             'bias_vs_z_std',
                                'growth function':        'growth_factor_de',
                                'pk-gauss-smoothing':     'pk_gauss_smoothing',
                                }

    lya_dla_cross_parameters_dict = {'bias_eta_DLA':              '1.     1.      None    None    fixed',
                                     'beta_DLA':                  '0.5    1.      None    None    free',
                                     'alpha_DLA':                 '0.0    0.      None    None    fixed',
                                     'par binsize LYA(LYA)xDLA':  '4      0.      None    None    fixed',
                                     'per binsize LYA(LYA)xDLA':  '4      0.      None    None    fixed',
                                     'drp_DLA':                   '0.     0.1     None    None    fixed',
                                     'sigma_velo_lorentz_DLA':    '0.     0.1     None    None    fixed',
                                     }

    lya_dla_cross_options_dict = {'data':         lya_dla_cross_data_dict,
                                  'cuts':         lya_dla_cross_cuts_dict,
                                  'model':        lya_dla_cross_model_dict,
                                  'parameters':   lya_dla_cross_parameters_dict,
                                  }

    config_filepath = fit_dir + '/config_lya_dla_cross_{}.ini'.format(suffix)
    make_config_file(config_filepath,lya_dla_cross_options_dict)
    configs += [config_filepath]

    chi2_filepath = fit_dir + 'chi2_{}_{}.ini'.format(name,suffix)
    result_filepath = fit_dir + 'result_{}_{}.h5'.format(name,suffix)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_lya_auto__lya_dla_cross__dla_auto_fit_files(fits_dir,exp_filepaths,rmin=20.,rmax=160.,afix='free'):

    name = 'lya_auto__lya_dla_cross__dla_auto'
    suffix = 'rmin{}_rmax{}_a{}'.format(rmin,rmax,afix)
    fit_dir = fits_dir + '/' + name + '/'
    submit_utils.check_dir(fit_dir)

    exp_files = [exp_filepaths[key] for key in exp_filepaths.keys() if key in name]
    zeff = get_zeff(exp_files)
    f = get_growth_rate(zeff,Om_z0=args.fid_Om)

    configs = []

    lya_auto_data_dict = {'name':            'LYA(LYA)xLYA(LYA)',
                          'tracer1':         'LYA',
                          'tracer2':         'LYA',
                          'tracer1-type':    'continuous',
                          'tracer2-type':    'continuous',
                          'filename':        exp_filepaths['lya_auto'],
                          'ell-max':         6,
                          }

    lya_auto_cuts_dict = {'rp-min':  0.,
                          'rp-max':  200.,
                          'rt-min':  0.,
                          'rt-max':  200.,
                          'r-min':   rmin,
                          'r-max':   rmax,
                          'mu-min':  0.,
                          'mu-max':  1.
                          }

    lya_auto_model_dict = {'model-pk':           'pk_kaiser',
                           'model-xi':           'xi',
                           'z evol LYA':         'bias_vs_z_std',
                           'growth function':    'growth_factor_de',
                           'pk-gauss-smoothing': 'pk_gauss_smoothing',
                           }

    lya_auto_parameters_dict = {'ap':                            '1.     0.1     0.8     1.2     {}'.format(afix),
                                'at':                            '1.     0.1     0.8     1.2     {}'.format(afix),
                                'bao_amp':                       '1.     0.      None    None    fixed',
                                'sigmaNL_per':                   '0.     0.      None    None    fixed',
                                'sigmaNL_par':                   '0.     0.      None    None    fixed',
                                'growth_rate':                   '{}     0.      None    None    free'.format(f),
                                'bias_eta_LYA':                  '-0.1   1.      None    None    free',
                                'beta_LYA':                      '0.5    1.      None    None    free',
                                'alpha_LYA':                     '2.9    0.      None    None    fixed',
                                'par binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'per binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'par_sigma_smooth':              '2.     2.      None    None    free',
                                'per_sigma_smooth':              '2.     2.      None    None    free',
                                }

    lya_auto_options_dict = {'data':         lya_auto_data_dict,
                             'cuts':         lya_auto_cuts_dict,
                             'model':        lya_auto_model_dict,
                             'parameters':   lya_auto_parameters_dict,
                             }

    config_filepath = fit_dir + '/config_lya_auto_{}.ini'.format(suffix)
    make_config_file(config_filepath,lya_auto_options_dict)
    configs += [config_filepath]

    lya_dla_cross_data_dict = {'name':            'LYA(LYA)xDLA',
                               'tracer1':         'LYA',
                               'tracer2':         'DLA',
                               'tracer1-type':    'continuous',
                               'tracer2-type':    'discrete',
                               'filename':        exp_filepaths['lya_dla_cross'],
                               'ell-max':         6,
                               }

    lya_dla_cross_cuts_dict = {'rp-min':  -200.,
                               'rp-max':  200.,
                               'rt-min':  0.,
                               'rt-max':  200.,
                               'r-min':   rmin,
                               'r-max':   rmax,
                               'mu-min':  -1.,
                               'mu-max':  1.
                               }

    lya_dla_cross_model_dict = {'model-pk':               'pk_kaiser',
                                'model-xi':               'xi_drp',
                                'z evol LYA':             'bias_vs_z_std',
                                'velocity dispersion':    'pk_velo_lorentz',
                                'z evol DLA':             'bias_vs_z_std',
                                'growth function':        'growth_factor_de',
                                'pk-gauss-smoothing':     'pk_gauss_smoothing',
                                }

    lya_dla_cross_parameters_dict = {'bias_eta_DLA':              '1.     1.      None    None    fixed',
                                     'beta_DLA':                  '0.5    1.      None    None    free',
                                     'alpha_DLA':                 '0.0    0.      None    None    fixed',
                                     'par binsize LYA(LYA)xDLA':  '4      0.      None    None    fixed',
                                     'per binsize LYA(LYA)xDLA':  '4      0.      None    None    fixed',
                                     'drp_DLA':                   '0.     0.1     None    None    fixed',
                                     'sigma_velo_lorentz_DLA':    '0.     0.1     None    None    fixed',
                                     }

    lya_dla_cross_options_dict = {'data':         lya_dla_cross_data_dict,
                                  'cuts':         lya_dla_cross_cuts_dict,
                                  'model':        lya_dla_cross_model_dict,
                                  'parameters':   lya_dla_cross_parameters_dict,
                                  }

    config_filepath = fit_dir + '/config_lya_dla_cross_{}.ini'.format(suffix)
    make_config_file(config_filepath,lya_dla_cross_options_dict)
    configs += [config_filepath]

    dla_auto_data_dict = {'name':            'DLAxDLA',
                          'tracer1':         'DLA',
                          'tracer2':         'DLA',
                          'tracer1-type':    'discrete',
                          'tracer2-type':    'discrete',
                          'filename':        exp_filepaths['dla_auto'],
                          'ell-max':         6,
                          }

    dla_auto_cuts_dict = {'rp-min':  0.,
                          'rp-max':  200.,
                          'rt-min':  0.,
                          'rt-max':  200.,
                          'r-min':   rmin,
                          'r-max':   rmax,
                          'mu-min':  0.,
                          'mu-max':  1.
                          }

    dla_auto_model_dict = {'model-pk':               'pk_kaiser',
                           'model-xi':               'xi_drp',
                           'velocity dispersion':    'pk_velo_lorentz',
                           'z evol DLA':             'bias_vs_z_std',
                           'growth function':        'growth_factor_de',
                           'pk-gauss-smoothing':     'pk_gauss_smoothing',
                           }

    dla_auto_parameters_dict = {'par binsize DLAxDLA':  '4      0.      None    None    fixed',
                                'per binsize DLAxDLA':  '4      0.      None    None    fixed',
                                }

    dla_auto_options_dict = {'data':         dla_auto_data_dict,
                             'cuts':         dla_auto_cuts_dict,
                             'model':        dla_auto_model_dict,
                             'parameters':   dla_auto_parameters_dict,
                             }

    config_filepath = fit_dir + '/config_dla_auto_{}.ini'.format(suffix)
    make_config_file(config_filepath,dla_auto_options_dict)
    configs += [config_filepath]

    chi2_filepath = fit_dir + 'chi2_{}_{}.ini'.format(name,suffix)
    result_filepath = fit_dir + 'result_{}_{}.h5'.format(name,suffix)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_qso_auto__qso_dla_cross__dla_auto_fit_files(fits_dir,exp_filepaths,rmin=20.,rmax=160.,afix='free'):

    name = 'qso_auto__qso_dla_cross__qso_auto'
    suffix = 'rmin{}_rmax{}_a{}'.format(rmin,rmax,afix)
    fit_dir = fits_dir + '/' + name + '/'
    submit_utils.check_dir(fit_dir)

    exp_files = [exp_filepaths[key] for key in exp_filepaths.keys() if key in name]
    zeff = get_zeff(exp_files)
    f = get_growth_rate(zeff,Om_z0=args.fid_Om)

    configs = []

    qso_auto_data_dict = {'name':            'QSOxQSO',
                          'tracer1':         'QSO',
                          'tracer2':         'QSO',
                          'tracer1-type':    'discrete',
                          'tracer2-type':    'discrete',
                          'filename':        exp_filepaths['qso_auto'],
                          'ell-max':         6,
                          }

    qso_auto_cuts_dict = {'rp-min':  0.,
                          'rp-max':  200.,
                          'rt-min':  0.,
                          'rt-max':  200.,
                          'r-min':   rmin,
                          'r-max':   rmax,
                          'mu-min':  0.,
                          'mu-max':  1.
                          }

    qso_auto_model_dict = {'model-pk':               'pk_kaiser',
                           'model-xi':               'xi_drp',
                           'velocity dispersion':    'pk_velo_lorentz',
                           'z evol QSO':             'bias_vs_z_std',
                           'growth function':        'growth_factor_de',
                           'pk-gauss-smoothing':     'pk_gauss_smoothing',
                           }

    qso_auto_parameters_dict = {'ap':                       '1.     0.1     0.8     1.2     {}'.format(afix),
                                'at':                       '1.     0.1     0.8     1.2     {}'.format(afix),
                                'bao_amp':                  '1.     0.      None    None    fixed',
                                'sigmaNL_per':              '0.     0.      None    None    fixed',
                                'sigmaNL_par':              '0.     0.      None    None    fixed',
                                'growth_rate':              '{}     0.      None    None    fixed'.format(f),
                                'bias_eta_QSO':             '1.0    1.      None    None    fixed',
                                'beta_QSO':                 '0.25   1.      None    None    free',
                                'alpha_QSO':                '1.44   0.      None    None    fixed',
                                'par binsize QSOxQSO':      '4      0.      None    None    fixed',
                                'per binsize QSOxQSO':      '4      0.      None    None    fixed',
                                'par_sigma_smooth':         '2.     2.      None    None    free',
                                'per_sigma_smooth':         '2.     2.      None    None    free',
                                'drp_QSO':                  '0.     0.1     None    None    fixed',
                                'sigma_velo_lorentz_QSO':   '0.     0.1     None    None    fixed',
                                }

    qso_auto_options_dict = {'data':         qso_auto_data_dict,
                             'cuts':         qso_auto_cuts_dict,
                             'model':        qso_auto_model_dict,
                             'parameters':   qso_auto_parameters_dict,
                             }

    config_filepath = fit_dir + '/config_qso_auto_{}.ini'.format(suffix)
    make_config_file(config_filepath,qso_auto_options_dict)
    configs += [config_filepath]

    qso_dla_cross_data_dict = {'name':            'QSOxDLA',
                               'tracer1':         'QSO',
                               'tracer2':         'DLA',
                               'tracer1-type':    'continuous',
                               'tracer2-type':    'discrete',
                               'filename':        exp_filepaths['lya_dla_cross'],
                               'ell-max':         6,
                               }

    qso_dla_cross_cuts_dict = {'rp-min':  -200.,
                               'rp-max':  200.,
                               'rt-min':  0.,
                               'rt-max':  200.,
                               'r-min':   rmin,
                               'r-max':   rmax,
                               'mu-min':  -1.,
                               'mu-max':  1.
                               }

    qso_dla_cross_model_dict = {'model-pk':               'pk_kaiser',
                                'model-xi':               'xi_drp',
                                'z evol LYA':             'bias_vs_z_std',
                                'velocity dispersion':    'pk_velo_lorentz',
                                'z evol DLA':             'bias_vs_z_std',
                                'growth function':        'growth_factor_de',
                                'pk-gauss-smoothing':     'pk_gauss_smoothing',
                                }

    qso_dla_cross_parameters_dict = {'bias_eta_DLA':              '1.     1.      None    None    fixed',
                                     'beta_DLA':                  '0.5    1.      None    None    free',
                                     'alpha_DLA':                 '0.0    0.      None    None    fixed',
                                     'par binsize QSOxDLA':       '4      0.      None    None    fixed',
                                     'per binsize QSOxDLA':       '4      0.      None    None    fixed',
                                     'drp_DLA':                   '0.     0.1     None    None    fixed',
                                     'sigma_velo_lorentz_DLA':    '0.     0.1     None    None    fixed',
                                     }

    qso_dla_cross_options_dict = {'data':         qso_dla_cross_data_dict,
                                  'cuts':         qso_dla_cross_cuts_dict,
                                  'model':        qso_dla_cross_model_dict,
                                  'parameters':   qso_dla_cross_parameters_dict,
                                  }

    config_filepath = fit_dir + '/config_qso_dla_cross_{}.ini'.format(suffix)
    make_config_file(config_filepath,qso_dla_cross_options_dict)
    configs += [config_filepath]

    dla_auto_data_dict = {'name':            'DLAxDLA',
                          'tracer1':         'DLA',
                          'tracer2':         'DLA',
                          'tracer1-type':    'discrete',
                          'tracer2-type':    'discrete',
                          'filename':        exp_filepaths['dla_auto'],
                          'ell-max':         6,
                          }

    dla_auto_cuts_dict = {'rp-min':  0.,
                          'rp-max':  200.,
                          'rt-min':  0.,
                          'rt-max':  200.,
                          'r-min':   rmin,
                          'r-max':   rmax,
                          'mu-min':  0.,
                          'mu-max':  1.
                          }

    dla_auto_model_dict = {'model-pk':               'pk_kaiser',
                           'model-xi':               'xi_drp',
                           'velocity dispersion':    'pk_velo_lorentz',
                           'z evol DLA':             'bias_vs_z_std',
                           'growth function':        'growth_factor_de',
                           'pk-gauss-smoothing':     'pk_gauss_smoothing',
                           }

    dla_auto_parameters_dict = {'par binsize DLAxDLA':  '4      0.      None    None    fixed',
                                'per binsize DLAxDLA':  '4      0.      None    None    fixed',
                                }

    dla_auto_options_dict = {'data':         dla_auto_data_dict,
                             'cuts':         dla_auto_cuts_dict,
                             'model':        dla_auto_model_dict,
                             'parameters':   dla_auto_parameters_dict,
                             }

    config_filepath = fit_dir + '/config_dla_auto_{}.ini'.format(suffix)
    make_config_file(config_filepath,dla_auto_options_dict)
    configs += [config_filepath]

    chi2_filepath = fit_dir + 'chi2_{}_{}.ini'.format(name,suffix)
    result_filepath = fit_dir + 'result_{}_{}.h5'.format(name,suffix)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_all_correlations_fit_files(fits_dir,exp_filepaths,rmin=20.,rmax=160.,afix='free'):

    name = 'all'
    suffix = 'rmin{}_rmax{}_a{}'.format(rmin,rmax,afix)
    fit_dir = fits_dir + '/' + name + '/'
    submit_utils.check_dir(fit_dir)

    exp_files = [exp_filepaths[key] for key in exp_filepaths.keys() if key in name]
    zeff = get_zeff(exp_files)
    f = get_growth_rate(zeff,Om_z0=args.fid_Om)

    configs = []

    lya_auto_data_dict = {'name':            'LYA(LYA)xLYA(LYA)',
                          'tracer1':         'LYA',
                          'tracer2':         'LYA',
                          'tracer1-type':    'continuous',
                          'tracer2-type':    'continuous',
                          'filename':        exp_filepaths['lya_auto'],
                          'ell-max':         6,
                          }

    lya_auto_cuts_dict = {'rp-min':  0.,
                          'rp-max':  200.,
                          'rt-min':  0.,
                          'rt-max':  200.,
                          'r-min':   rmin,
                          'r-max':   rmax,
                          'mu-min':  0.,
                          'mu-max':  1.
                          }

    lya_auto_model_dict = {'model-pk':           'pk_kaiser',
                           'model-xi':           'xi',
                           'z evol LYA':         'bias_vs_z_std',
                           'growth function':    'growth_factor_de',
                           'pk-gauss-smoothing': 'pk_gauss_smoothing',
                           }

    lya_auto_parameters_dict = {'ap':                            '1.     0.1     0.8     1.2     {}'.format(afix),
                                'at':                            '1.     0.1     0.8     1.2     {}'.format(afix),
                                'bao_amp':                       '1.     0.      None    None    fixed',
                                'sigmaNL_per':                   '0.     0.      None    None    fixed',
                                'sigmaNL_par':                   '0.     0.      None    None    fixed',
                                'growth_rate':                   '{}     0.      None    None    free'.format(f),
                                'bias_eta_LYA':                  '-0.1   1.      None    None    free',
                                'beta_LYA':                      '0.5    1.      None    None    free',
                                'alpha_LYA':                     '2.9    0.      None    None    fixed',
                                'par binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'per binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'par_sigma_smooth':              '2.     2.      None    None    free',
                                'per_sigma_smooth':              '2.     2.      None    None    free',
                                }

    lya_auto_options_dict = {'data':         lya_auto_data_dict,
                             'cuts':         lya_auto_cuts_dict,
                             'model':        lya_auto_model_dict,
                             'parameters':   lya_auto_parameters_dict,
                             }

    config_filepath = fit_dir + '/config_lya_auto_{}.ini'.format(suffix)
    make_config_file(config_filepath,lya_auto_options_dict)
    configs += [config_filepath]

    lya_qso_cross_data_dict = {'name':            'LYA(LYA)xQSO',
                               'tracer1':         'LYA',
                               'tracer2':         'QSO',
                               'tracer1-type':    'continuous',
                               'tracer2-type':    'discrete',
                               'filename':        exp_filepaths['lya_qso_cross'],
                               'ell-max':         6,
                               }

    lya_qso_cross_cuts_dict = {'rp-min':  -200.,
                               'rp-max':  200.,
                               'rt-min':  0.,
                               'rt-max':  200.,
                               'r-min':   rmin,
                               'r-max':   rmax,
                               'mu-min':  -1.,
                               'mu-max':  1.
                               }

    lya_qso_cross_model_dict = {'model-pk':               'pk_kaiser',
                                'model-xi':               'xi_drp',
                                'z evol LYA':             'bias_vs_z_std',
                                'velocity dispersion':    'pk_velo_lorentz',
                                'z evol QSO':             'bias_vs_z_std',
                                'growth function':        'growth_factor_de',
                                'pk-gauss-smoothing':     'pk_gauss_smoothing',
                                }

    lya_qso_cross_parameters_dict = {'bias_eta_QSO':              '1.     1.      None    None    fixed',
                                     'beta_QSO':                  '0.25   1.      None    None    free',
                                     'alpha_QSO':                 '1.44   0.      None    None    fixed',
                                     'par binsize LYA(LYA)xQSO':  '4      0.      None    None    fixed',
                                     'per binsize LYA(LYA)xQSO':  '4      0.      None    None    fixed',
                                     'drp_QSO':                   '0.     0.1     None    None    fixed',
                                     'sigma_velo_lorentz_QSO':    '0.     0.1     None    None    fixed',
                                     }

    lya_qso_cross_options_dict = {'data':         lya_qso_cross_data_dict,
                                  'cuts':         lya_qso_cross_cuts_dict,
                                  'model':        lya_qso_cross_model_dict,
                                  'parameters':   lya_qso_cross_parameters_dict,
                                  }

    config_filepath = fit_dir + '/config_lya_qso_cross_{}.ini'.format(suffix)
    make_config_file(config_filepath,lya_qso_cross_options_dict)
    configs += [config_filepath]

    qso_auto_data_dict = {'name':            'QSOxQSO',
                          'tracer1':         'QSO',
                          'tracer2':         'QSO',
                          'tracer1-type':    'discrete',
                          'tracer2-type':    'discrete',
                          'filename':        exp_filepaths['qso_auto'],
                          'ell-max':         6,
                          }

    qso_auto_cuts_dict = {'rp-min':  0.,
                          'rp-max':  200.,
                          'rt-min':  0.,
                          'rt-max':  200.,
                          'r-min':   rmin,
                          'r-max':   rmax,
                          'mu-min':  0.,
                          'mu-max':  1.
                          }

    qso_auto_model_dict = {'model-pk':               'pk_kaiser',
                           'model-xi':               'xi_drp',
                           'velocity dispersion':    'pk_velo_lorentz',
                           'z evol QSO':             'bias_vs_z_std',
                           'growth function':        'growth_factor_de',
                           'pk-gauss-smoothing':     'pk_gauss_smoothing',
                           }

    qso_auto_parameters_dict = {'par binsize QSOxQSO':  '4      0.      None    None    fixed',
                                'per binsize QSOxQSO':  '4      0.      None    None    fixed',
                                }

    qso_auto_options_dict = {'data':         qso_auto_data_dict,
                             'cuts':         qso_auto_cuts_dict,
                             'model':        qso_auto_model_dict,
                             'parameters':   qso_auto_parameters_dict,
                             }

    config_filepath = fit_dir + '/config_qso_auto_{}.ini'.format(suffix)
    make_config_file(config_filepath,qso_auto_options_dict)
    configs += [config_filepath]

    qso_dla_cross_data_dict = {'name':            'QSOxDLA',
                               'tracer1':         'QSO',
                               'tracer2':         'DLA',
                               'tracer1-type':    'continuous',
                               'tracer2-type':    'discrete',
                               'filename':        exp_filepaths['lya_dla_cross'],
                               'ell-max':         6,
                               }

    qso_dla_cross_cuts_dict = {'rp-min':  -200.,
                               'rp-max':  200.,
                               'rt-min':  0.,
                               'rt-max':  200.,
                               'r-min':   rmin,
                               'r-max':   rmax,
                               'mu-min':  -1.,
                               'mu-max':  1.
                               }

    qso_dla_cross_model_dict = {'model-pk':               'pk_kaiser',
                                'model-xi':               'xi_drp',
                                'z evol LYA':             'bias_vs_z_std',
                                'velocity dispersion':    'pk_velo_lorentz',
                                'z evol DLA':             'bias_vs_z_std',
                                'growth function':        'growth_factor_de',
                                'pk-gauss-smoothing':     'pk_gauss_smoothing',
                                }

    qso_dla_cross_parameters_dict = {'bias_eta_DLA':              '1.     1.      None    None    fixed',
                                     'beta_DLA':                  '0.5    1.      None    None    free',
                                     'alpha_DLA':                 '0.0    0.      None    None    fixed',
                                     'par binsize QSOxDLA':       '4      0.      None    None    fixed',
                                     'per binsize QSOxDLA':       '4      0.      None    None    fixed',
                                     'drp_DLA':                   '0.     0.1     None    None    fixed',
                                     'sigma_velo_lorentz_DLA':    '0.     0.1     None    None    fixed',
                                     }

    qso_dla_cross_options_dict = {'data':         qso_dla_cross_data_dict,
                                  'cuts':         qso_dla_cross_cuts_dict,
                                  'model':        qso_dla_cross_model_dict,
                                  'parameters':   qso_dla_cross_parameters_dict,
                                  }

    config_filepath = fit_dir + '/config_qso_dla_cross_{}.ini'.format(suffix)
    make_config_file(config_filepath,qso_dla_cross_options_dict)
    configs += [config_filepath]

    dla_auto_data_dict = {'name':            'DLAxDLA',
                          'tracer1':         'DLA',
                          'tracer2':         'DLA',
                          'tracer1-type':    'discrete',
                          'tracer2-type':    'discrete',
                          'filename':        exp_filepaths['dla_auto'],
                          'ell-max':         6,
                          }

    dla_auto_cuts_dict = {'rp-min':  0.,
                          'rp-max':  200.,
                          'rt-min':  0.,
                          'rt-max':  200.,
                          'r-min':   rmin,
                          'r-max':   rmax,
                          'mu-min':  0.,
                          'mu-max':  1.
                          }

    dla_auto_model_dict = {'model-pk':               'pk_kaiser',
                           'model-xi':               'xi_drp',
                           'velocity dispersion':    'pk_velo_lorentz',
                           'z evol DLA':             'bias_vs_z_std',
                           'growth function':        'growth_factor_de',
                           'pk-gauss-smoothing':     'pk_gauss_smoothing',
                           }

    dla_auto_parameters_dict = {'par binsize DLAxDLA':  '4      0.      None    None    fixed',
                                'per binsize DLAxDLA':  '4      0.      None    None    fixed',
                                }

    dla_auto_options_dict = {'data':         dla_auto_data_dict,
                             'cuts':         dla_auto_cuts_dict,
                             'model':        dla_auto_model_dict,
                             'parameters':   dla_auto_parameters_dict,
                             }

    config_filepath = fit_dir + '/config_dla_auto_{}.ini'.format(suffix)
    make_config_file(config_filepath,dla_auto_options_dict)
    configs += [config_filepath]

    lya_dla_cross_data_dict = {'name':            'LYA(LYA)xDLA',
                               'tracer1':         'LYA',
                               'tracer2':         'DLA',
                               'tracer1-type':    'continuous',
                               'tracer2-type':    'discrete',
                               'filename':        exp_filepaths['lya_dla_cross'],
                               'ell-max':         6,
                               }

    lya_dla_cross_cuts_dict = {'rp-min':  -200.,
                               'rp-max':  200.,
                               'rt-min':  0.,
                               'rt-max':  200.,
                               'r-min':   rmin,
                               'r-max':   rmax,
                               'mu-min':  -1.,
                               'mu-max':  1.
                               }

    lya_dla_cross_model_dict = {'model-pk':               'pk_kaiser',
                                'model-xi':               'xi_drp',
                                'z evol LYA':             'bias_vs_z_std',
                                'velocity dispersion':    'pk_velo_lorentz',
                                'z evol DLA':             'bias_vs_z_std',
                                'growth function':        'growth_factor_de',
                                'pk-gauss-smoothing':     'pk_gauss_smoothing',
                                }

    lya_dla_cross_parameters_dict = {'par binsize LYA(LYA)xDLA':  '4      0.      None    None    fixed',
                                     'per binsize LYA(LYA)xDLA':  '4      0.      None    None    fixed',
                                     }

    lya_dla_cross_options_dict = {'data':         lya_dla_cross_data_dict,
                                  'cuts':         lya_dla_cross_cuts_dict,
                                  'model':        lya_dla_cross_model_dict,
                                  'parameters':   lya_dla_cross_parameters_dict,
                                  }

    config_filepath = fit_dir + '/config_lya_dla_cross_{}.ini'.format(suffix)
    make_config_file(config_filepath,lya_dla_cross_options_dict)
    configs += [config_filepath]

    chi2_filepath = fit_dir + 'chi2_{}_{}.ini'.format(name,suffix)
    result_filepath = fit_dir + 'result_{}_{}.h5'.format(name,suffix)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

################################################################################

#Calculate the effective redshift.
def get_zeff(fi,rmin=80.,rmax=120.):
    zeffs = []
    weights = []
    for f in fi:
        h = fits.open(f)
        R = np.sqrt(h[1].data['RP']**2 + h[1].data['RT']**2)
        cells = (R > rmin) * (R < rmax)
        #Should this be weighted by WE or by NB?
        zeff = np.average(h[1].data['Z'][cells],weights=h[1].data['NB'][cells])
        weight = np.sum(h[1].data['NB'][cells])
        h.close()
        zeffs += [zeff]
        weights += [weight]
    zeff = np.average(zeffs,weights=weights)
    return zeff

#Calculate f. Assume Or=0 for now.
# TODO: include the input Om and Or from picca
def get_growth_rate(z,Om_z0=0.3147):
    Om = Om_z0 * ((1+z)**3) / (Om_z0 * ((1+z)**3) + 1 - Om_z0)
    Ol = 1 - Om
    f = (Om**0.6) + (Ol/70.)*(1 + Om/2.)
    return f

#Calculate beta_QSO by interpolating the input bias.
# TODO: work out how to set the location more appropriately.
def get_beta_obj(z,obj='QSO'):
    b_qso_data_loc='/global/homes/j/jfarr/Projects/LyaCoLoRe/input_files/Bz_qso_G18.txt'
    b_qso_data = np.loadtxt(b_qso_data_loc)
    if obj == 'QSO':
        bias_of_z = interp1d(b_qso_data[:,0],b_qso_data[:,1])
    elif obj == 'DLA':
        print(' -> -> -> DLA bias assumed to be constant with z.')
        bias_of_z = interp1d(b_qso_data[:,0],2.*np.ones_like(b_qso_data[:,0]))
    bias = bias_of_z(z)
    f = get_growth_rate(z)
    beta = f/bias
    print(z,bias,f,beta)
    return beta

################################################################################

a_dir = args.base_dir+'/analysis/'
submit_utils.check_dir(a_dir)
ac_dir = a_dir+'/correlation_functions/'
submit_utils.check_dir(ac_dir)

vers = []
for v_rea in args.v_realisations:
    ver = 'v{}.{}.{}'.format(args.v_maj,args.v_min,v_rea)
    vers += [ver]
if args.fit_stack:
    vers += ['stack']

for ver in vers:
    print('\nMaking fit files for version {}:'.format(ver))

    acv_dir = ac_dir+'/'+ver+'/'
    submit_utils.check_dir(acv_dir)
    acvf_dir = acv_dir+'/fits/'
    submit_utils.check_dir(acvf_dir)
    acvm_dir = acv_dir+'/measurements/'
    submit_utils.check_dir(acvm_dir)

    exp_filepaths = {'lya_auto':        acvm_dir + '/lya_auto/correlations/cf_exp_lya_auto.fits.gz',
                     'qso_auto':        acvm_dir + '/qso_auto/correlations/co_exp_qso_auto.fits.gz',
                     'dla_auto':        acvm_dir + '/dla_auto/correlations/co_exp_dla_auto.fits.gz',
                     'lya_aa_auto':     acvm_dir + '/lya_aa_auto/correlations/cf_exp_lya_aa_auto.fits.gz',
                     'lya_qso_cross':   acvm_dir + '/lya_qso_cross/correlations/xcf_exp_lya_qso_cross.fits.gz',
                     'lya_dla_cross':   acvm_dir + '/lya_dla_cross/correlations/xcf_exp_lya_dla_cross.fits.gz',
                     #'qso_dla_cross':   acvm_dir + '/qso_dla_cross/correlations/co_exp_lya_dla_cross.fits.gz',
                     }

    #For each correlation, make the correlation information file
    for name,exp_filepath in exp_filepaths.items():
        if name == 'lya_auto':
            corr_type = 'cf'
            q1 = q2 = 'LYA'
        if name == 'qso_auto':
            corr_type = 'co'
            q1 = q2 = 'QSO'
        if name == 'dla_auto':
            corr_type = 'co'
            q1 = q2 = 'DLA'
        if name == 'lya_aa_auto':
            corr_type = 'cf'
            q1 = q2 = 'LYA'
        if name == 'lya_qso_cross':
            corr_type = 'xcf'
            q1 = 'LYA'
            q2 = 'QSO'
        if name == 'lya_dla_cross':
            corr_type = 'xcf'
            q1 = 'LYA'
            q2 = 'DLA'
        if name == 'qso_dla_cross':
            corr_type = 'co'
            q1 = 'QSO'
            q2 = 'DLA'
        make_corr_info_file(acvm_dir+'/'+name+'/corr_info.txt',corr_type,q1,q2)

    #create config files for doing fits
    for rmin in args.rmin_values:
        for rmax in args.rmax_values:
            for afix in args.afix_values:
                print(' -> making files for rmin={}, rmax={}, afix={}'.format(rmin,rmax,afix))

                if args.fit_lya_auto:
                    print(' -> -> making lya_auto files')
                    make_lya_auto_fit_files(acvf_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_qso_auto:
                    print(' -> -> making qso_auto files')
                    make_qso_auto_fit_files(acvf_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_dla_auto:
                    print(' -> -> making dla_auto files')
                    make_dla_auto_fit_files(acvf_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_aa_auto:
                    print(' -> -> making lya_aa_auto files')
                    make_lya_aa_auto_fit_files(acvf_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_qso_cross:
                    print(' -> -> making lya_qso_cross files')
                    make_lya_qso_cross_fit_files(acvf_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_dla_cross:
                    print(' -> -> making lya_dla_cross files')
                    make_lya_dla_cross_fit_files(acvf_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_qso_dla_cross:
                    print(' -> -> making qso_dla_cross files')
                    make_qso_dla_cross_fit_files(acvf_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_auto__lya_qso_cross:
                    print(' -> -> making joint lya_auto__lya_qso_cross files')
                    make_lya_auto__lya_qso_cross_fit_files(acvf_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_auto__lya_qso_cross__qso_auto:
                    print(' -> -> making joint lya_auto__lya_qso_cross__qso_auto files')
                    make_lya_auto__lya_qso_cross__qso_auto_fit_files(acvf_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_auto__lya_dla_cross:
                    print(' -> -> making joint lya_auto__lya_dla_cross files')
                    make_lya_auto__lya_dla_cross_fit_files(acvf_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_auto__lya_dla_cross__dla_auto:
                    print(' -> -> making joint lya_auto__lya_dla_cross__dla_auto files')
                    make_lya_auto__lya_dla_cross__dla_auto_fit_files(acvf_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_qso_auto__qso_dla_cross__dla_auto:
                    print(' -> -> making joint qso_auto__qso_dla_cross__dla_auto files')
                    make_qso_auto__qso_dla_cross__dla_auto_fit_files(acvf_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_all_correlations:
                    print(' -> -> making joint all correlations files')
                    make_all_correlations_fit_files(acvf_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
