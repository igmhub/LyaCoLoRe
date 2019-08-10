#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#locations and names

parser.add_argument('--base-dir', type = str, default = None, required=True,
                    help = 'base directory')

parser.add_argument('--out-dir', type = str, default = None, required=False,
                    help = 'output directory')

parser.add_argument('--parameter-filename', type = str, default = 'parameters.txt', required=False,
                    help = 'name of parameter file')

#cf descriptors

parser.add_argument('--corr-type', type = str, default = None, required=True,
                    help = 'base directory')

parser.add_argument('--quantity', type = str, default = None, required=True,
                    help = 'base directory')

parser.add_argument('--npixels', type = int, default = None, required=True,
                    help = 'number of pixels')

parser.add_argument('--quant-code', type = str, default = None, required=True,
                    help = 'base directory')

parser.add_argument('--zmin', type = float, default = None, required=True,
                    help = 'minimum redshift cut value')

parser.add_argument('--zmax', type = float, default = None, required=True,
                    help = 'maximum redshift cut value')

parser.add_argument('--corr-filename', type = str, default = 'cf.fits.gz', required=False,
                    help = 'name of correlation file')

parser.add_argument('--corr-exp-filename', type = str, default = 'cf_exp.fits.gz', required=False,
                    help = 'name of exported correlation file')

parser.add_argument('--nside', type = int, default = 16, required=False,
                    help = 'HEALPix nside for output files (must be 2^n)')

parser.add_argument('--rpmin', type = float, default = 0., required=False,
                    help = 'minimum parallel separation')

parser.add_argument('--rpmax', type = float, default = 160., required=False,
                    help = 'maximum parallel separation')

parser.add_argument('--rtmin', type = float, default = 0., required=False,
                    help = 'minimum transverse separation')

parser.add_argument('--rtmax', type = float, default = 160., required=False,
                    help = 'maximum transverse separation')

parser.add_argument('--np', type = int, default = 40, required=False,
                    help = 'number of parallel bins')

parser.add_argument('--nt', type = int, default = 40, required=False,
                    help = 'number of transverse bins')

parser.add_argument('--include-metals', action="store_true", default = False, required=False,
                    help = 'include metal absorption')

#Fit variable

parser.add_argument('--rmin-values', type = float, default = [0.], required=True,
                    help = 'minimum separation', nargs='*')

parser.add_argument('--rmax-values', type = float, default = [160.], required=True,
                    help = 'maximum separation', nargs='*')

parser.add_argument('--afix-values', type = str, default = None, required=True,
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

parser.add_argument('--fit-lya-auto--lya-qso-cross--qso--auto', action="store_true", default = False, required=False,
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
    chi2_text += '[output]\n'
    chi2_text += 'filename = {}\n\n'.format(result_filename)
    chi2_text += '[fiducial]\n'
    chi2_text += 'filename = /PlanckDR12/PlanckDR12.fits\n'
    file = open(filepath,'w')
    file.write(chi2_text)
    file.close()

    return

def make_parameter_file(filepath,args):

    parameter_file_text = ''
    parameter_file_text += 'correl_type = {}\n'.format(args.corr_type)
    parameter_file_text += 'quantity = {}\n'.format(args.quantity)
    parameter_file_text += 'N_side = {}\n'.format(args.nside)
    parameter_file_text += 'N_pixels = {}\n'.format(args.npixels)
    parameter_file_text += 'quantities = {}\n'.format(args.quant_code)
    parameter_file_text += 'rpmin = {}\n'.format(args.rpmin)
    parameter_file_text += 'rpmax = {}\n'.format(args.rpmax)
    parameter_file_text += 'rtmin = {}\n'.format(args.rtmin)
    parameter_file_text += 'rtmax = {}\n'.format(args.rtmax)
    parameter_file_text += 'np = {}\n'.format(args.np)
    parameter_file_text += 'nt = {}\n'.format(args.nt)
    parameter_file_text += 'zmin = {}\n'.format(args.zmin)
    parameter_file_text += 'zmax = {}\n'.format(args.zmax)
    parameter_file_text += 'zeff = {}\n'.format(zeff)
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
        print(' ')

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

    parameters_dict = {'ap':                            '1.     0.1     0.5     1.5     {}'.format(afix),
                       'at':                            '1.     0.1     0.5     1.5     {}'.format(afix),
                       'bao_amp':                       '1.     0.      None    None    fixed',
                       'sigmaNL_per':                   '0.     0.      None    None    fixed',
                       'sigmaNL_par':                   '0.     0.      None    None    fixed',
                       'growth_rate':                   '{}     0.      None    None    fixed'.format(f), # TODO: NEED TO calculate
                       'bias_eta_LYA':                  '-0.1   1.      None    None    free',
                       'beta_LYA':                      '0.5    1.      None    None    free',
                       'alpha_LYA':                     '2.9    0.      None    None    fixed',
                       'par binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                       'per binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                       'par_sigma_smooth':              '1.42   0.4     None    None    fixed',
                       'per_sigma_smooth':              '3.78   0.4     None    None    fixed',
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

    parameters_dict = {'ap':                        '1.     0.1     0.5     1.5     {}'.format(afix),
                       'at':                        '1.     0.1     0.5     1.5     {}'.format(afix),
                       'bao_amp':                   '1.     0.      None    None    fixed',
                       'sigmaNL_per':               '0.     0.      None    None    fixed',
                       'sigmaNL_par':               '0.     0.      None    None    fixed',
                       'growth_rate':               '{}     0.      None    None    fixed'.format(), # TODO: NEED TO calculate
                       'bias_eta_QSO':              '1.     1.      None    None    fixed',
                       'beta_QSO':                  '0.25   1.      None    None    free',
                       'alpha_QSO':                 '1.44   0.      None    None    fixed',
                       'par binsize QSOxQSO':       '4      0.      None    None    fixed',
                       'per binsize QSOxQSO':       '4      0.      None    None    fixed',
                       'par_sigma_smooth':          '1.42   0.4     None    None    fixed',
                       'per_sigma_smooth':          '3.78   0.4     None    None    fixed',
                       'drp_QSO':                   '0.     0.1     None    None    fixed',
                       'sigma_velo_lorentz_QSO':    '2.     0.1     None    None    free',
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

    parameters_dict = {'ap':                        '1.     0.1     0.5     1.5     {}'.format(afix),
                       'at':                        '1.     0.1     0.5     1.5     {}'.format(afix),
                       'bao_amp':                   '1.     0.      None    None    fixed',
                       'sigmaNL_per':               '0.     0.      None    None    fixed',
                       'sigmaNL_par':               '0.     0.      None    None    fixed',
                       'growth_rate':               '{}     0.      None    None    fixed'.format(), # TODO: NEED TO calculate
                       'bias_eta_DLA':              '1.     1.      None    None    fixed',
                       'beta_DLA':                  '0.5    1.      None    None    free',
                       'alpha_DLA':                 '0.0    0.      None    None    fixed',
                       'par binsize DLAxDLA':       '4      0.      None    None    fixed',
                       'per binsize DLAxDLA':       '4      0.      None    None    fixed',
                       'par_sigma_smooth':          '1.42   0.4     None    None    fixed',
                       'per_sigma_smooth':          '3.78   0.4     None    None    fixed',
                       'drp_DLA':                   '0.     0.1     None    None    fixed',
                       'sigma_velo_lorentz_DLA':    '2.     0.1     None    None    free',
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

    parameters_dict = {'ap':                            '1.     0.1     0.5     1.5     {}'.format(afix),
                       'at':                            '1.     0.1     0.5     1.5     {}'.format(afix),
                       'bao_amp':                       '1.     0.      None    None    fixed',
                       'sigmaNL_per':                   '0.     0.      None    None    fixed',
                       'sigmaNL_par':                   '0.     0.      None    None    fixed',
                       'growth_rate':                   '{}     0.      None    None    fixed'.format(), # TODO: NEED TO calculate
                       'bias_eta_LYA':                  '-0.1   1.      None    None    free',
                       'beta_LYA':                      '0.5    1.      None    None    free',
                       'alpha_LYA':                     '2.9    0.      None    None    fixed',
                       'par binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                       'per binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                       'par_sigma_smooth':              '1.42   0.4     None    None    fixed',
                       'per_sigma_smooth':              '3.78   0.4     None    None    fixed',
                       }

    metals_dict = {'filename':      '{}/data/additional_data/metal_dmat_z_0_10.fits\n'.format(args.base_dir),
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

    parameters_dict = {'ap':                        '1.     0.1     0.5     1.5     {}'.format(afix),
                       'at':                        '1.     0.1     0.5     1.5     {}'.format(afix),
                       'bao_amp':                   '1.     0.      None    None    fixed',
                       'sigmaNL_per':               '0.     0.      None    None    fixed',
                       'sigmaNL_par':               '0.     0.      None    None    fixed',
                       'growth_rate':               '{}     0.      None    None    fixed'.format(), # TODO: NEED TO calculate
                       'bias_eta_LYA':              '-0.1   1.      None    None    free',
                       'beta_LYA':                  '0.5    1.      None    None    free',
                       'alpha_LYA':                 '2.9    0.      None    None    fixed',
                       'bias_eta_QSO':              '1.     1.      None    None    fixed',
                       'beta_QSO':                  '{}     1.      None    None    fixed'.format(), # TODO: NEED TO calculate
                       'alpha_QSO':                 '1.44   0.      None    None    fixed',
                       'par binsize LYA(LYA)xQSO':  '4      0.      None    None    fixed',
                       'per binsize LYA(LYA)xQSO':  '4      0.      None    None    fixed',
                       'par_sigma_smooth':          '1.42   0.4     None    None    fixed',
                       'per_sigma_smooth':          '3.78   0.4     None    None    fixed',
                       'drp_QSO':                   '0.     0.1     None    None    fixed',
                       'sigma_velo_lorentz_QSO':    '2.     0.1     None    None    free',
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

    parameters_dict = {'ap':                        '1.     0.1     0.5     1.5     {}'.format(afix),
                       'at':                        '1.     0.1     0.5     1.5     {}'.format(afix),
                       'bao_amp':                   '1.     0.      None    None    fixed',
                       'sigmaNL_per':               '0.     0.      None    None    fixed',
                       'sigmaNL_par':               '0.     0.      None    None    fixed',
                       'growth_rate':               '{}     0.      None    None    fixed'.format(), # TODO: NEED TO calculate
                       'bias_eta_LYA':              '-0.1   1.      None    None    free',
                       'beta_LYA':                  '0.5    1.      None    None    free',
                       'alpha_LYA':                 '2.9    0.      None    None    fixed',
                       'bias_eta_DLA':              '1.     1.      None    None    fixed',
                       'beta_DLA':                  '{}     1.      None    None    fixed'.format(), # TODO: NEED TO calculate
                       'alpha_DLA':                 '0.0    0.      None    None    fixed',
                       'par binsize LYA(LYA)xDLA':  '4      0.      None    None    fixed',
                       'per binsize LYA(LYA)xDLA':  '4      0.      None    None    fixed',
                       'par_sigma_smooth':          '1.42   0.4     None    None    fixed',
                       'per_sigma_smooth':          '3.78   0.4     None    None    fixed',
                       'drp_DLA':                   '0.     0.1     None    None    fixed',
                       'sigma_velo_lorentz_DLA':    '2.     0.1     None    None    free',
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

    parameters_dict = {'ap':                        '1.     0.1     0.5     1.5     {}'.format(afix),
                       'at':                        '1.     0.1     0.5     1.5     {}'.format(afix),
                       'bao_amp':                   '1.     0.      None    None    fixed',
                       'sigmaNL_per':               '0.     0.      None    None    fixed',
                       'sigmaNL_par':               '0.     0.      None    None    fixed',
                       'growth_rate':               '{}     0.      None    None    fixed'.format(), # TODO: NEED TO calculate
                       'bias_eta_QSO':              '1.     1.      None    None    fixed',
                       'beta_QSO':                  '0.5    1.      None    None    free', # TODO: can we leave both betas free?
                       'alpha_QSO':                 '1.44   0.      None    None    fixed',
                       'bias_eta_DLA':              '1.     1.      None    None    fixed',
                       'beta_DLA':                  '0.25   1.      None    None    free', # TODO: can we leave both betas free?
                       'alpha_DLA':                 '0.0    0.      None    None    fixed',
                       'par binsize QSOxDLA':       '4      0.      None    None    fixed',
                       'per binsize QSOxDLA':       '4      0.      None    None    fixed',
                       'par_sigma_smooth':          '1.42   0.4     None    None    fixed',
                       'per_sigma_smooth':          '3.78   0.4     None    None    fixed',
                       'drp_QSO':                   '0.     0.1     None    None    fixed',
                       'sigma_velo_lorentz_QSO':    '2.     0.1     None    None    free',
                       'drp_DLA':                   '0.     0.1     None    None    fixed',
                       'sigma_velo_lorentz_DLA':    '2.     0.1     None    None    free',
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

    lya_auto_parameters_dict = {'ap':                            '1.     0.1     0.5     1.5     {}'.format(afix),
                                'at':                            '1.     0.1     0.5     1.5     {}'.format(afix),
                                'bao_amp':                       '1.     0.      None    None    fixed',
                                'sigmaNL_per':                   '0.     0.      None    None    fixed',
                                'sigmaNL_par':                   '0.     0.      None    None    fixed',
                                'growth_rate':                   '{}     0.      None    None    fixed'.format(), # TODO: NEED TO calculate
                                'bias_eta_LYA':                  '-0.1   1.      None    None    free',
                                'beta_LYA':                      '0.5    1.      None    None    free',
                                'alpha_LYA':                     '2.9    0.      None    None    fixed',
                                'par binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'per binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'par_sigma_smooth':              '1.42   0.4     None    None    fixed',
                                'per_sigma_smooth':              '3.78   0.4     None    None    fixed',
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
                                     'beta_QSO':                  '{}     1.      None    None    fixed'.format(), # TODO: NEED TO calculate
                                     'alpha_QSO':                 '1.44   0.      None    None    fixed',
                                     'par binsize LYA(LYA)xQSO':  '4      0.      None    None    fixed',
                                     'per binsize LYA(LYA)xQSO':  '4      0.      None    None    fixed',
                                     'drp_QSO':                   '0.     0.1     None    None    fixed',
                                     'sigma_velo_lorentz_QSO':    '2.     0.1     None    None    free',
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

    lya_auto_parameters_dict = {'ap':                            '1.     0.1     0.5     1.5     {}'.format(afix),
                                'at':                            '1.     0.1     0.5     1.5     {}'.format(afix),
                                'bao_amp':                       '1.     0.      None    None    fixed',
                                'sigmaNL_per':                   '0.     0.      None    None    fixed',
                                'sigmaNL_par':                   '0.     0.      None    None    fixed',
                                'growth_rate':                   '{}     0.      None    None    fixed'.format(), # TODO: NEED TO calculate
                                'bias_eta_LYA':                  '-0.1   1.      None    None    free',
                                'beta_LYA':                      '0.5    1.      None    None    free',
                                'alpha_LYA':                     '2.9    0.      None    None    fixed',
                                'par binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'per binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'par_sigma_smooth':              '1.42   0.4     None    None    fixed',
                                'per_sigma_smooth':              '3.78   0.4     None    None    fixed',
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
                                     'sigma_velo_lorentz_QSO':    '2.     0.1     None    None    free',
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

    lya_auto_parameters_dict = {'ap':                            '1.     0.1     0.5     1.5     {}'.format(afix),
                                'at':                            '1.     0.1     0.5     1.5     {}'.format(afix),
                                'bao_amp':                       '1.     0.      None    None    fixed',
                                'sigmaNL_per':                   '0.     0.      None    None    fixed',
                                'sigmaNL_par':                   '0.     0.      None    None    fixed',
                                'growth_rate':                   '{}     0.      None    None    fixed'.format(), # TODO: NEED TO calculate
                                'bias_eta_LYA':                  '-0.1   1.      None    None    free',
                                'beta_LYA':                      '0.5    1.      None    None    free',
                                'alpha_LYA':                     '2.9    0.      None    None    fixed',
                                'par binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'per binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'par_sigma_smooth':              '1.42   0.4     None    None    fixed',
                                'per_sigma_smooth':              '3.78   0.4     None    None    fixed',
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
                                     'beta_DLA':                  '{}     1.      None    None    fixed'.format(), # TODO: NEED TO calculate
                                     'alpha_DLA':                 '1.44   0.      None    None    fixed',
                                     'par binsize LYA(LYA)xDLA':  '4      0.      None    None    fixed',
                                     'per binsize LYA(LYA)xDLA':  '4      0.      None    None    fixed',
                                     'drp_DLA':                   '0.     0.1     None    None    fixed',
                                     'sigma_velo_lorentz_DLA':    '2.     0.1     None    None    free',
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

    lya_auto_parameters_dict = {'ap':                            '1.     0.1     0.5     1.5     {}'.format(afix),
                                'at':                            '1.     0.1     0.5     1.5     {}'.format(afix),
                                'bao_amp':                       '1.     0.      None    None    fixed',
                                'sigmaNL_per':                   '0.     0.      None    None    fixed',
                                'sigmaNL_par':                   '0.     0.      None    None    fixed',
                                'growth_rate':                   '{}     0.      None    None    fixed'.format(), # TODO: NEED TO calculate
                                'bias_eta_LYA':                  '-0.1   1.      None    None    free',
                                'beta_LYA':                      '0.5    1.      None    None    free',
                                'alpha_LYA':                     '2.9    0.      None    None    fixed',
                                'par binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'per binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'par_sigma_smooth':              '1.42   0.4     None    None    fixed',
                                'per_sigma_smooth':              '3.78   0.4     None    None    fixed',
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
                                     'sigma_velo_lorentz_DLA':    '2.     0.1     None    None    free',
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

    qso_auto_parameters_dict = {'ap':                   '1.     0.1     0.5     1.5     {}'.format(afix),
                                'at':                   '1.     0.1     0.5     1.5     {}'.format(afix),
                                'bao_amp':              '1.     0.      None    None    fixed',
                                'sigmaNL_per':          '0.     0.      None    None    fixed',
                                'sigmaNL_par':          '0.     0.      None    None    fixed',
                                'growth_rate':          '{}     0.      None    None    fixed'.format(), # TODO: NEED TO calculate
                                'bias_eta_QSO':         '1.0    1.      None    None    fixed',
                                'beta_QSO':             '0.25   1.      None    None    free',
                                'alpha_QSO':            '1.44   0.      None    None    fixed',
                                'par binsize QSOxQSO':  '4      0.      None    None    fixed',
                                'per binsize QSOxQSO':  '4      0.      None    None    fixed',
                                'par_sigma_smooth':     '1.42   0.4     None    None    fixed',
                                'per_sigma_smooth':     '3.78   0.4     None    None    fixed',
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
                                     'sigma_velo_lorentz_DLA':    '2.     0.1     None    None    free',
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

    lya_auto_parameters_dict = {'ap':                            '1.     0.1     0.5     1.5     {}'.format(afix),
                                'at':                            '1.     0.1     0.5     1.5     {}'.format(afix),
                                'bao_amp':                       '1.     0.      None    None    fixed',
                                'sigmaNL_per':                   '0.     0.      None    None    fixed',
                                'sigmaNL_par':                   '0.     0.      None    None    fixed',
                                'growth_rate':                   '{}     0.      None    None    fixed'.format(), # TODO: NEED TO calculate
                                'bias_eta_LYA':                  '-0.1   1.      None    None    free',
                                'beta_LYA':                      '0.5    1.      None    None    free',
                                'alpha_LYA':                     '2.9    0.      None    None    fixed',
                                'par binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'per binsize LYA(LYA)xLYA(LYA)': '4      0.      None    None    fixed',
                                'par_sigma_smooth':              '1.42   0.4     None    None    fixed',
                                'per_sigma_smooth':              '3.78   0.4     None    None    fixed',
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
                                     'sigma_velo_lorentz_QSO':    '2.     0.1     None    None    free',
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
                                     'sigma_velo_lorentz_DLA':    '2.     0.1     None    None    free',
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
def get_beta_qso(z,b_qso_of_z_loc='/global/homes/j/jfarr/Projects/LyaCoLoRe/input_files/Bz_qso_G18.txt'):
    b_qso_of_z = np.loadtxt(b_qso_of_z_loc)
    b_qso = np.interp(z,b_qso_of_z[:,0],b_qso_of_z[:,1])
    f = get_growth_rate(z)
    beta_qso = f/b_qso
    return beta_qso

################################################################################

exp_filepaths = {'lya_auto':        acvm_dir + '/lya_auto/cf_lya_auto_exp.fits.gz',
                 'qso_auto':        acvm_dir + '/qso_auto/co_qso_auto_exp.fits.gz',
                 'dla_auto':        acvm_dir + '/lya_auto/co_dla_auto_exp.fits.gz',
                 'lya_aa_auto':     acvm_dir + '/lya_aa_auto/cf_lya_aa_auto_exp.fits.gz',
                 'lya_qso_cross':   acvm_dir + '/lya_qso_cross/xcf_lya_qso_cross_exp.fits.gz',
                 'lya_dla_cross':   acvm_dir + '/lya_dla_cross/xcf_lya_dla_cross_exp.fits.gz',
                 'qso_dla_cross':   acvm_dir + '/qso_dla_cross/co_lya_dla_cross_exp.fits.gz',
                 }

#For each correlation, make the parameter file
#for name,exp_filepath in exp_filepaths.items():
#    make_parameter_file(acvm_dir+'/'+name+'/parameters.txt',args)

for v_rea in args.v_realisations:

    ac_dir = a_dir+'/correlation_functions/'
    submit_utils.check_dir(ac_dir)
    acv_dir = ac_dir+'/'+ver+'/'
    submit_utils.check_dir(acv_dir)
    acvf_dir = acv_dir+'/fits/'
    submit_utils.check_dir(acvf_dir)

    #create config files for doing fits
    for rmin in args.rmin_values:
        for rmax in args.rmax_values:
            for afix in args.afix_values:

                if args.fit_lya_auto:
                    make_lya_auto_fit_files(fits_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_qso_auto:
                    make_qso_auto_fit_files(fits_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_dla_auto:
                    make_dla_auto_fit_files(fits_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_aa_auto:
                    make_lya_aa_auto_fit_files(fits_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_qso_cross:
                    make_lya_qso_cross_fit_files(fits_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_dla_cross:
                    make_lya_dla_cross_fit_files(fits_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_qso_dla_cross:
                    make_qso_dla_cross_fit_files(fits_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_auto__lya_qso_cross:
                    make_lya_auto__lya_qso_cross_fit_files(fits_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_auto__lya_qso_cross__qso_auto:
                    make_lya_auto__lya_qso_cross__qso_auto_fit_files(fits_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_auto__lya_dla_cross:
                    make_lya_auto__lya_dla_cross_fit_files(fits_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_auto__lya_dla_cross__dla_auto:
                    make_lya_auto__lya_dla_cross__dla_auto_fit_files(fits_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_qso_auto__qso_dla_cross__dla_auto:
                    make_qso_auto__qso_dla_cross__dla_auto_fit_files(fits_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_all_correlations:
                    make_all_correlations_fit_files(fits_dir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
