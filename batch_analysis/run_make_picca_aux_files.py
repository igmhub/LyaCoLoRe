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

args = parser.parse_args()

################################################################################

if args.out_dir is None:
    args.out_dir = args.base_dir

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
def make_lya_auto_fit_files(exp_filepaths,config_filepaths,parameter_filepath,chi2_filepath,result_filepath,rmin=20.,rmax=160.,afix='free'):

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
                       'growth_rate':                   '{}     0.      None    None    fixed'.format(), # TODO: NEED TO calculate
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

    make_config_file(config_filepaths['lya_auto'],options_dict)
    configs += [config_filepaths['lya_auto']]

    make_parameter_filepath(parameter_filepath,args)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return options_dict

def make_qso_auto_fit_files(exp_filepaths,config_filepaths,parameter_filepath,chi2_filepath,result_filepath,rmin=20.,rmax=160.,afix='free'):

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

    make_config_file(config_filepaths['qso_auto'],options_dict)
    configs += [config_filepaths['qso_auto']]

    make_parameter_filepath(parameter_filepath,args)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_dla_auto_fit_files(exp_filepaths,config_filepaths,parameter_filepath,chi2_filepath,result_filepath,rmin=20.,rmax=160.,afix='free'):

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

    make_config_file(config_filepaths['dla_auto'],options_dict)
    configs += [config_filepaths['dla_auto']]

    make_parameter_filepath(parameter_filepath,args)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_lya_aa_auto_fit_files(exp_filepaths,config_filepaths,parameter_filepath,chi2_filepath,result_filepath,rmin=20.,rmax=160.,afix='free'):

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

    make_config_file(config_filepaths['lya_aa_auto'],options_dict)
    configs += [config_filepaths['lya_aa_auto']]

    make_parameter_filepath(parameter_filepath,args)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_lya_qso_cross_fit_files(exp_filepaths,config_filepaths,parameter_filepath,chi2_filepath,result_filepath,rmin=20.,rmax=160.,afix='free'):

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

    make_config_file(config_filepaths['lya_qso_cross'],options_dict)
    configs += [config_filepaths['lya_qso_cross']]

    make_parameter_filepath(parameter_filepath,args)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_lya_dla_cross_fit_files(exp_filepaths,config_filepaths,parameter_filepath,chi2_filepath,result_filepath,rmin=20.,rmax=160.,afix='free'):

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

    make_config_file(config_filepaths['lya_dla_cross'],options_dict)
    configs += [config_filepaths['lya_dla_cross']]

    make_parameter_filepath(parameter_filepath,args)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_qso_dla_cross_fit_files(exp_filepaths,config_filepaths,parameter_filepath,chi2_filepath,result_filepath,rmin=20.,rmax=160.,afix='free'):

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

    make_config_file(config_filepaths['qso_dla_cross'],options_dict)
    configs += [config_filepaths['qso_dla_cross']]

    make_parameter_filepath(parameter_filepath,args)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_lya_auto__lya_qso_cross_fit_files(exp_filepaths,config_filepaths,parameter_filepath,chi2_filepath,result_filepath,rmin=20.,rmax=160.,afix='free'):

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

    make_config_file(config_filepaths['lya_auto'],lya_auto_options_dict)
    configs += [config_filepaths['lya_auto']]

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

    make_config_file(config_filepaths['lya_qso_cross'],lya_qso_cross_options_dict)
    configs += [config_filepaths['lya_qso_cross']]

    # TODO: What to do with this
    make_parameter_filepath(parameter_filepath,args)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_lya_auto__lya_qso_cross__qso_auto_fit_files(exp_filepaths,config_filepaths,parameter_filepath,chi2_filepath,result_filepath,rmin=20.,rmax=160.,afix='free'):

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

    make_config_file(config_filepaths['lya_auto'],lya_auto_options_dict)
    configs += [config_filepaths['lya_auto']]

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

    make_config_file(config_filepaths['lya_qso_cross'],lya_qso_cross_options_dict)
    configs += [config_filepaths['lya_qso_cross']]

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

    make_config_file(config_filepaths['lya_qso_cross'],lya_qso_cross_options_dict)
    configs += [config_filepaths['lya_qso_cross']]

    # TODO: What to do with this
    make_parameter_filepath(parameter_filepath,args)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_lya_auto__lya_dla_cross_fit_files(exp_filepaths,config_filepaths,parameter_filepath,chi2_filepath,result_filepath,rmin=20.,rmax=160.,afix='free'):

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

    make_config_file(config_filepaths['lya_auto'],lya_auto_options_dict)
    configs += [config_filepaths['lya_auto']]

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

    make_config_file(config_filepaths['lya_dla_cross'],lya_dla_cross_options_dict)
    configs += [config_filepaths['lya_dla_cross']]

    # TODO: What to do with this
    make_parameter_filepath(parameter_filepath,args)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_lya_auto__lya_dla_cross__dla_auto_fit_files(exp_filepaths,config_filepaths,parameter_filepath,chi2_filepath,result_filepath,rmin=20.,rmax=160.,afix='free'):

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

    make_config_file(config_filepaths['lya_auto'],lya_auto_options_dict)
    configs += [config_filepaths['lya_auto']]

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

    make_config_file(config_filepaths['lya_dla_cross'],lya_dla_cross_options_dict)
    configs += [config_filepaths['lya_dla_cross']]

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

    make_config_file(config_filepaths['dla_auto'],dla_auto_options_dict)
    configs += [config_filepaths['dla_auto']]

    # TODO: What to do with this
    make_parameter_filepath(parameter_filepath,args)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_qso_auto__qso_dla_cross__dla_auto_fit_files(exp_filepaths,config_filepaths,parameter_filepath,chi2_filepath,result_filepath,rmin=20.,rmax=160.,afix='free'):

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

    make_config_file(config_filepaths['qso_auto'],qso_auto_options_dict)
    configs += [config_filepaths['qso_auto']]

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

    make_config_file(config_filepaths['qso_dla_cross'],qso_dla_cross_options_dict)
    configs += [config_filepaths['qso_dla_cross']]

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

    make_config_file(config_filepaths['dla_auto'],dla_auto_options_dict)
    configs += [config_filepaths['dla_auto']]

    # TODO: What to do with this
    make_parameter_filepath(parameter_filepath,args)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

def make_all_correlations_fit_files(exp_filepaths,config_filepaths,parameter_filepath,chi2_filepath,result_filepath,rmin=20.,rmax=160.,afix='free'):

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

    make_config_file(config_filepaths['lya_auto'],lya_auto_options_dict)
    configs += [config_filepaths['lya_auto']]

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

    make_config_file(config_filepaths['lya_qso_cross'],lya_qso_cross_options_dict)
    configs += [config_filepaths['lya_qso_cross']]

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

    make_config_file(config_filepaths['lya_qso_cross'],lya_qso_cross_options_dict)
    configs += [config_filepaths['lya_qso_cross']]





    # TODO: What to do with this
    make_parameter_filepath(parameter_filepath,args)
    make_chi2_file(chi2_filepath,zeff,configs,result_filepath)

    return

################################################################################

#Calculate the effective redshift.
h = fits.open(args.base_dir+'/'+args.corr_exp_filename)
R = np.sqrt(h[1].data['RP']**2 + h[1].data['RT']**2)
cells = (R > 80.) * (R < 120.)
zeff = np.average(h[1].data['Z'][cells],weights=h[1].data['NB'][cells])
h.close()

#Calculate f. Assume Or=0 for now.
# TODO: include the input Om and Or from picca.
Om_z0 = 0.3147
Om = Om_z0 * ((1+zeff)**3) / (Om_z0 * ((1+zeff)**3) + 1 - Om_z0)
Ol = 1 - Om
f = (Om**0.6) + (Ol/70.)*(1 + Om/2.)

#Calculate beta_QSO by interpolating the input bias.
b_qso_of_z = np.loadtxt('/global/homes/j/jfarr/Projects/LyaCoLoRe/input_files/Bz_qso_G18.txt')
b_qso = np.interp(zeff,b_qso_of_z[:,0],b_qso_of_z[:,1])
beta_qso = f/b_qso

#Calculate beta_DLA by interpolating the input bias.
# TODO: Want to read b_dla from the input ideally.
b_dla = 2.
beta_dla = f/b_qso

#Write parameter file.
make_parameter_file(args.base_dir+args.parameter_filename,args)

#create config files for doing fits
for rmin in args.rmin_values:
    for rmax in args.rmax_values:
        for afix in args.afix_values:

            suffix = '{}r_a{}'.format(int(rmin),afix)
            config_filename = args.base_dir+'config_{}_{}.ini'.format(args.corr_type,suffix)
            result_filename = args.base_dir+'/result_{}.h5\n\n'.format(suffix)
            chi2_filename = args.base_dir+'chi2_{}.ini'.format(suffix)

            make_chi2_file(chi2_filename,zeff,[config_filename],result_filename)

            config_text = ''
            config_text += '#!/bin/bash -l\n\n'
            config_text += '[data]\n'
            if args.corr_type == 'cf':
                config_text += 'name = LYA(LYA)xLYA(LYA)\n'
                config_text += 'tracer1 = LYA\n'
                config_text += 'tracer2 = LYA\n'
                config_text += 'tracer1-type = continuous\n'
                config_text += 'tracer2-type = continuous\n'
            elif args.corr_type == 'xcf':
                config_text += 'name = LYA(LYA)xQSO\n'
                config_text += 'tracer1 = QSO\n'
                config_text += 'tracer2 = LYA\n'
                config_text += 'tracer1-type = discrete\n'
                config_text += 'tracer2-type = continuous\n'
            elif args.corr_type == 'co':
                config_text += 'name = QSOxQSO\n'
                config_text += 'tracer1 = QSO\n'
                config_text += 'tracer2 = QSO\n'
                config_text += 'tracer1-type = discrete\n'
                config_text += 'tracer2-type = discrete\n'
            config_text += 'filename = {}\n'.format(args.base_dir+args.corr_exp_filename)
            config_text += 'ell-max = 6\n\n'
            config_text += '[cuts]\n'
            config_text += 'rp-min = {}\n'.format(args.rpmin)
            config_text += 'rp-max = {}\n\n'.format(args.rpmax)
            config_text += 'rt-min = {}\n'.format(args.rtmin)
            config_text += 'rt-max = {}\n\n'.format(args.rtmax)
            config_text += 'r-min = {}\n'.format(rmin)
            config_text += 'r-max = {}\n\n'.format(rmax)
            if args.corr_type == 'cf':
                config_text += 'mu-min = 0.\n'
                config_text += 'mu-max = +1.\n\n'
            elif args.corr_type == 'xcf':
                config_text += 'mu-min = -1.\n'
                config_text += 'mu-max = +1.\n\n'
            elif args.corr_type == 'co':
                config_text += 'mu-min = -1.\n'
                config_text += 'mu-max = +1.\n\n'
            config_text += '[model]\n'
            config_text += 'model-pk = pk_kaiser\n'
            if args.corr_type == 'cf':
                config_text += 'model-xi = xi\n'
                config_text += 'z evol LYA = bias_vs_z_std\n'
            elif args.corr_type == 'xcf':
                config_text += 'model-xi = xi_drp\n'
                config_text += 'z evol LYA = bias_vs_z_std\n'
                config_text += 'z evol QSO = bias_vs_z_std\n'
                config_text += 'velocity dispersion = pk_velo_lorentz\n'
            elif args.corr_type == 'co':
                config_text += 'model-xi = xi_drp\n'
                config_text += 'z evol QSO = bias_vs_z_std\n'
                config_text += 'velocity dispersion = pk_velo_lorentz\n'
            config_text += 'growth function = growth_factor_de\n'
            config_text += 'pk-gauss-smoothing = pk_gauss_smoothing\n\n'
            config_text += '[parameters]\n\n'
            if args.corr_type == 'cf':
                config_text += 'drp_QSO                = 0. 0.1   None None fixed\n'
            elif args.corr_type == 'xcf':
                config_text += 'drp_QSO                = 0. 0.1   None None free\n'
            elif args.corr_type == 'co':
                config_text += 'drp_QSO                = 0. 0.1   None None fixed\n'
            if args.corr_type == 'xcf' or args.corr_type == 'co':
                config_text += 'sigma_velo_lorentz_QSO = 2. 0.1    None None free\n\n'
            config_text += 'ap = 1. 0.1 0.5 1.5 {}\n'.format(afix)
            config_text += 'at = 1. 0.1 0.5 1.5 {}\n'.format(afix)
            config_text += 'bao_amp = 1. 0. None None fixed\n\n'
            config_text += 'sigmaNL_per = 0.     0. None None fixed\n'
            config_text += 'sigmaNL_par = 0.     0.1 None None fixed\n'
            config_text += 'growth_rate = {} 0. None None fixed\n\n'.format(f)
            if args.corr_type == 'cf':
                config_text += 'bias_eta_LYA  = -0.0003512  1. None None free\n'
                config_text += 'beta_LYA  = 0.5    1. None None free\n'
                config_text += 'alpha_LYA = 2.9     0. None None fixed\n\n'
            if args.corr_type == 'xcf':
                config_text += 'bias_eta_LYA  = -0.0003512  1. None None free\n'
                config_text += 'beta_LYA  = 1.4    1. None None free\n'
                config_text += 'alpha_LYA = 2.9     0. None None fixed\n'
                config_text += 'bias_eta_QSO  = 1. 0. None None fixed\n'
                config_text += 'beta_QSO      = {} 0.1 None None fixed\n'.format(beta_QSO)
                config_text += 'alpha_QSO = 1.44     0. None None fixed\n\n'
            if args.corr_type == 'co':
                config_text += 'bias_eta_QSO  = 1. 1. None None fixed\n'
                config_text += 'beta_QSO      = {} 1. None None free\n'.format(beta_QSO)
                config_text += 'alpha_QSO = 1.44     0. None None fixed\n\n'
            if args.corr_type == 'cf':
                config_text += 'par binsize LYA(LYA)xLYA(LYA) = 4 0. None None fixed\n'
                config_text += 'per binsize LYA(LYA)xLYA(LYA) = 4 0. None None fixed\n\n'
            elif args.corr_type == 'xcf':
                config_text += 'par binsize LYA(LYA)xQSO = 4 0. None None fixed\n'
                config_text += 'per binsize LYA(LYA)xQSO = 4 0. None None fixed\n\n'
            elif args.corr_type == 'co':
                config_text += 'par binsize QSOxQSO = 4 0. None None fixed\n'
                config_text += 'per binsize QSOxQSO = 4 0. None None fixed\n\n'
            config_text += 'par_sigma_smooth = 1.42 0.4 None None fixed\n'
            config_text += 'per_sigma_smooth = 3.78 0.4 None None fixed\n\n'
            if args.include_metals:
                config_text += '[metals]\n'
                config_text += 'filename = {}/metal_dmat_z_0_10.fits\n'.format(args.base_dir)
                config_text += 'model-pk-met = pk_kaiser\n'
                config_text += 'model-xi-met = cached_xi_kaiser\n'
                config_text += 'z evol = bias_vs_z_std\n'
                config_text += 'in tracer1 = SiII(1260) SiIII(1207) SiII(1193) SiII(1190)\n'
                config_text += 'in tracer2 = SiII(1260) SiIII(1207) SiII(1193) SiII(1190)\n'
            file = open(config_filename,'w')
            file.write(config_text)
            file.close()
