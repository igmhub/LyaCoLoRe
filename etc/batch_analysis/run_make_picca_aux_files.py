#!/usr/bin/env python

import argparse
from astropy.io import fits
import numpy as np
import os
from scipy.interpolate import interp1d

from lyacolore import submit_utils

################################################################################

tracer_types = {
    'lya': 'continuous',
    'qso': 'discrete',
    'dla': 'discrete',
}

beta_qso_fix_dict = {
    'lyalya_qso': 'fixed',
    'lyalyb_qso': 'fixed',
    'qso_qso': 'free',
    'qso_dla': 'free',
    'lyalya_lyalya__lyalya_qso': 'free',
    'lyalya_lyalyb__lyalyb_qso': 'free',
    'lyalya_qso__lyalyb_qso': 'fixed',
    'lyalya_lyalya__lyalya_lyalyb__lyalya_qso__lyalyb_qso': 'free',
}
beta_dla_fix_dict = {
    'lyalya_dla': 'fixed',
    'lyalyb_dla': 'fixed',
    'qso_dla': 'fixed',
    'dla_dla': 'free',
    'lyalya_lyalya__lyalya_dla': 'free',
    'lyalya_lyalyb__lyalyb_dla': 'free',
    'lyalya_dla__lyalyb_dla': 'fixed',
    'lyalya_lyalya__lyalya_lyalyb__lyalya_dla__lyalyb_dla': 'free',
}

fits_recognised = [
    'lyalya_lyalya',
    'lyalya_lyalyb',
    'lyalya_qso',
    'lyalyb_qso',
    'lyalya_dla',
    'lyalyb_dla',
    'qso_qso',
    'qso_dla',
    'dla_dla',
    'lyalya_lyalya__lyalya_lyalyb',
    'lyalya_lyalya__lyalya_qso',
    'lyalya_lyalyb__lyalyb_qso',
    'lyalya_lyalya__lyalya_dla',
    'lyalya_lyalyb__lyalyb_dla',
    'lyalya_qso__lyalyb_qso',
    'lyalya_dla__lyalyb_dla',
    'lyalya_lyalya__lyalya_lyalyb__lyalya_qso__lyalyb_qso',
    'lyalya_lyalya__lyalya_lyalyb__lyalya_dla__lyalyb_dla',
]

################################################################################

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#locations and names

parser.add_argument('--corr-dirs', type = str, nargs='*', default = None, required=True,
                    help = 'base directory')

parser.add_argument('--fid-Om', type=float, default=0.315, required=False,
                    help='Omega_matter(z=0) of fiducial LambdaCDM cosmology')

parser.add_argument('--fid-Or', type=float, default=0., required=False,
                    help='Omega_radiation(z=0) of fiducial LambdaCDM cosmology')

#Fit variables

parser.add_argument('--rmin-values', type = float, default = [40.], required=False,
                    help = 'minimum separation', nargs='*')

parser.add_argument('--rmax-values', type = float, default = [160.], required=False,
                    help = 'maximum separation', nargs='*')

parser.add_argument('--afix-values', type = str, default = ['free'], required=False,
                    help = 'maximum separation', nargs='*')

#Correlations to fit

parser.add_argument('--fitnames', type = str, default = ['lyalya_lyalya'], required=False,
                    help = 'fits to make, formatted as {corr1}__{corr2}__{corr3}\
                     where each of the "corr"s is a correlation name such as\
                     lyalya_lyalya or lyalya_qso', nargs='*')

args = parser.parse_args()

for f in args.fitnames:
    if f not in fits_recognised:
        raise ValueError('Fit {} not recognised!\nPlease choose from these options:\n{}'.format(f,fits_recognised))

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

## Function to make chi2 and config files for any given fit.
def make_fit_files(fitsdir,fitname,correlations,rmin=20.,rmax=180.,afix='free',beta_qso_fix='free',beta_dla_fix='free'):

    suffix = 'rmin{}_rmax{}_a{}'.format(rmin,rmax,afix)
    fitdir = os.path.join(fitsdir,fitname)
    submit_utils.check_dir(fitdir)

    exp_files = [c['exp_filepath'] for c in correlations.values()]
    zeff = get_zeff(exp_files)
    f = get_growth_rate(zeff,Om_z0=args.fid_Om)

    info = {cname: {'filepath': '', 'options_dict': {}} for cname in correlations.keys}

    ## For each correlation in the fit:
    configs = []
    config_filepaths = []
    for k,c in correlations.items():

        ## Make the data dict
        info[k]['options_dict']['data'] = {
            'name': k,
            'tracer1': c['tracer1'],
            'tracer1-type': c['tracer1-type'],
            'tracer2': c['tracer2'],
            'tracer2-type': c['tracer2-type'],
            'filename': c['exp_filepath'],
            'ell-max': 6,
        }

        ## Make the cuts dict
        info[k]['options_dict']['cuts'] = {
            'rp-min':  -200.,
            'rp-max':  200.,
            'rt-min':  0.,
            'rt-max':  200.,
            'r-min':   rmin,
            'r-max':   rmax,
            'mu-min':  -1.,
            'mu-max':  1.
        }

        ## Make the model dict
        info[k]['options_dict']['model'] = {
            'model-pk': 'pk_kaiser',
            'growth function': 'growth_factor_de',
            'pk-gauss-smoothing': 'pk_gauss_smoothing'
        }
        for t in set([c['tracer1'], c['tracer2']]):
            info[k]['options_dict']['model']['z evol {}'.format(t)] = 'bias_vs_z_std'
        if 'discrete' in set([c['tracer1-type'], c['tracer2-type']]):
            info[k]['options_dict']['model']['model-xi'] = 'xi_drp'
            info[k]['options_dict']['model']['velocity dispersion'] = 'pk_velo_lorentz'
        else:
            info[k]['options_dict']['model']['model-xi'] = 'xi'

        ## Make the parameters dict
        info[k]['options_dict']['parameters'] = {
            'ap':                       '1.     0.1     0.8     1.2     {}'.format(afix),
            'at':                       '1.     0.1     0.8     1.2     {}'.format(afix),
            'bao_amp':                  '1.     0.      None    None    fixed',
            'sigmaNL_per':              '0.     0.      None    None    fixed',
            'sigmaNL_par':              '0.     0.      None    None    fixed',
            'growth_rate':              '{}     0.      None    None    fixed'.format(f),
            'par binsize {}'.format(k): '4      0.      None    None    fixed',
            'per binsize {}'.format(k): '4      0.      None    None    fixed',
            'par_sigma_smooth':         '2.     2.      None    None    free',
            'per_sigma_smooth':         '2.     2.      None    None    free',
        }
        for t in set([c['tracer1'], c['tracer2']]):
            info[k]['options_dict']['parameters'] = {
                **info[k]['options_dict']['parameters'],
                **get_tracer_parameter_dict(t,beta_qso_fix=beta_qso_fix,beta_dla_fix=beta_dla_fix),
            }

        ## Make sure that parameters are not repeated from previous config files
        for k_past in configs:
            for p in info[k]['options_dict']['parameters'].keys():
                if p in info[k_past]['options_dict']['parameters'].keys():
                    _ = info[k]['options_dict']['parameters'].pop(p)

        ## Make the config file
        config_filepath = os.path.join(fitdir,'config_{}_{}.ini'.format(k,suffix))
        make_config_file(config_filepath,info[k]['options_dict'])
        configs += [k]
        config_filepaths += [config_filepath]

    ## Make the chi2 file
    chi2_filepath = os.path.join(fitdir,'chi2_{}_{}.ini'.format(fitname,suffix))
    result_filepath = os.path.join(fitdir,'result_{}_{}.h5'.format(fitname,suffix))
    make_chi2_file(chi2_filepath,zeff,config_filepaths,result_filepath)

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
    b_qso_data_loc='/global/homes/j/jfarr/Projects/LyaCoLoRe/etc/run_CoLoRe/input_files/Bz_qso_G18.txt'
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

def get_tracer_parameter_dict(tracer,beta_qso_fix='free',beta_dla_fix='free'):

    tracer_parameters = {
        'lya': {
            'bias_eta_lya':  '-0.2   0.1     None    None    free',
            'beta_lya':      '1.5    0.1     None    None    free',
            'alpha_lya':     '2.9    0.1     None    None    fixed',
        },
        'qso': {
            'bias_eta_qso':             '1.0    0.1     None    None    fixed',
            'beta_qso':                 '0.25   0.1     None    None    {}'.format(beta_qso_fix),
            'alpha_qso':                '1.44   0.1     None    None    fixed',
            'drp_qso':                  '0.     0.1     None    None    fixed',
            'sigma_velo_lorentz_qso':   '0.     0.1     None    None    fixed',
        },
        'dla': {
            'bias_eta_dla':             '1.0    0.1     None    None    fixed',
            'beta_dla':                 '0.25   0.1     None    None    {}'.format(beta_dla_fix),
            'alpha_dla':                '0.     0.1     None    None    fixed',
            'drp_dla':                  '0.     0.1     None    None    fixed',
            'sigma_velo_lorentz_dla':   '0.     0.1     None    None    fixed',
        }
    }

    return tracer_parameters[tracer]

################################################################################

for corr_dir in args.corr_dirs:
    print('\nMaking fit files for correlations in {}:'.format(corr_dir))

    dirloc = os.path.dirname(corr_dir)
    dirname = os.path.basename(corr_dir)
    analysis_dir = submit_utils.AnalysisDir(dirloc, dirname)

    correlations = {'lyalya_lyalya':    {'type': 'cf'},
                    'lyalya_lyalyb':    {'type': 'cf'},
                    'lyalya_qso':       {'type': 'xcf'},
                    'lyalyb_qso':       {'type': 'xcf'},
                    'lyalya_dla':       {'type': 'xcf'},
                    'lyalyb_dla':       {'type': 'xcf'},
                    'qso_qso':          {'type': 'co'},
                    'dla_dla':          {'type': 'co'},
                    'qso_dla':          {'type': 'co'},
                    }

    for k,c in correlations.items():
        for i,t in enumerate(k.split('_')):
            if t in ['lyalya','lyalyb']:
                c['tracer{}'.format(i)] = 'lya'
            else:
                c['tracer{}'.format(i)] = t
        c['tracer1-type'] = tracer_types[c['tracer1']]
        c['tracer2-type'] = tracer_types[c['tracer2']]
        c['exp_filepath'] = os.path.join(analysis_dir.corrdir,
                                        name,
                                        'measurements',
                                        '{}_exp_{}.fits.gz'.format(c['type'],k)
                                        )

    #create config files for doing fits
    for rmin in args.rmin_values:
        for rmax in args.rmax_values:
            for afix in args.afix_values:
                print(' -> making files for rmin={}, rmax={}, afix={}'.format(rmin,rmax,afix))
                for fitname in args.fitnames:

                    if 'qso' in fitname:
                        beta_qso_fix = beta_qso_fix_dict[fitname]
                    else:
                        beta_qso_fix = 'fixed'

                    if 'dla' in fitname:
                        beta_dla_fix = beta_dla_fix_dict[fitname]
                    else:
                        beta_dla_fix = 'fixed'

                    fit_correlations = {correlations[c] for c in fitname.split('__')}
                    print(' -> -> making {} files'.format(fitname))
                    make_fit_files(analysis_dir.fitsdir,fitname,fit_correlations,rmin=rmin,rmax=rmax,afix=afix,beta_qso_fix=beta_qso_fix,beta_dla_fix=beta_dla_fix)



"""



parser.add_argument('--fit-lyalya-lyalya', action="store_true", default = False, required=False,
                    help = 'make files for fitting the auto correlation of lya absorption in lya region')

parser.add_argument('--fit-lyalya-lyalyb', action="store_true", default = False, required=False,
                    help = 'make files for fitting the cross correlation of lya absorption in lya and lyb regions')

parser.add_argument('--fit-qso-qso', action="store_true", default = False, required=False,
                    help = 'make files for fitting the auto correlation of qsos')

parser.add_argument('--fit-dla-dla', action="store_true", default = False, required=False,
                    help = 'make files for fitting the auto correlation of dlas')

parser.add_argument('--fit-lyalya-qso', action="store_true", default = False, required=False,
                    help = 'make files for fitting the cross correlation of lya absorption in lya region and qsos')

parser.add_argument('--fit-lyalyb-qso', action="store_true", default = False, required=False,
                    help = 'make files for fitting the cross correlation of lya absorption in lyb region and qsos')

parser.add_argument('--fit-lyalya-dla', action="store_true", default = False, required=False,
                    help = 'make files for fitting the cross correlation of lya absorption in lya region and dlas')

parser.add_argument('--fit-lyalyb-dla', action="store_true", default = False, required=False,
                    help = 'make files for fitting the cross correlation of lya absorption in lyb region and dlas')

parser.add_argument('--fit-qso-dla', action="store_true", default = False, required=False,
                    help = 'make files for fitting the cross correlation of qsos and dlas')

parser.add_argument('--fit-full-auto', action="store_true", default = False, required=False,
                    help = 'make files for jointly fitting lyalya-lyalya and lyalya-lyalyb (DR16 "full auto")')

parser.add_argument('--fit-full-cross', action="store_true", default = False, required=False,
                    help = 'make files for jointly fitting lyalya-qso and lyalyb-qso (DR16 "full cross")')

parser.add_argument('--fit-all-combined', action="store_true", default = False, required=False,
                    help = 'make files for jointly fitting lyalya-lyalya, lyalya-lyalyb, lyalya-qso, and lyalyb-qso (DR16 "all combined")')

parser.add_argument('--fit-lyalya--lyalya-qso', action="store_true", default = False, required=False,
                    help = 'make files for joint lyalya-lyalya, lyalya-qso fit')

parser.add_argument('--fit-lyalya--lyalya-qso--qso-qso', action="store_true", default = False, required=False,
                    help = 'make files for joint lyalya-lyalya, lyalya-qso, qso-qso fit')

parser.add_argument('--fit-lyalya--lyalya-dla', action="store_true", default = False, required=False,
                    help = 'make files for joint lyalya-lyalya, lyalya-dla fit')

parser.add_argument('--fit-lyalya--lyalya-dla--dla-dla', action="store_true", default = False, required=False,
                    help = 'make files for joint lyalya-lyalya, lyalya-dla, dla-dla fit')

parser.add_argument('--fit-qso-qso--qso-dla--dla-dla', action="store_true", default = False, required=False,
                    help = 'make files for joint qso-qso, qso-dla, dla-dla fit')

parser.add_argument('--fit-all', action="store_true", default = False, required=False,
                    help = 'make files for fitting all correlations')

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

    if rmin >= 21:
        fix_SiIII1207_bias_eta = 'fixed'
        print('SiIII(1207) bias fixed to DR12 value as rmin is too large to fit it')
    else:
        fix_SiIII1207_bias_eta = 'free'
    if rmin >= 56:
        fix_SiII1193_bias_eta = 'fixed'
        print('SiII(1193) bias fixed to DR12 value as rmin is too large to fit it')
    else:
        fix_SiII1193_bias_eta = 'free'
    if rmin >= 64:
        fix_SiII1190_bias_eta = 'fixed'
        print('SiII(1190) bias fixed to DR12 value as rmin is too large to fit it')
    else:
        fix_SiII1190_bias_eta = 'free'
    if rmin >= 111:
        fix_SiII1260_bias_eta = 'fixed'
        print('SiII(1260) bias fixed to DR12 value as rmin is too large to fit it')
    else:
        fix_SiII1260_bias_eta = 'free'

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
                       'bias_eta_SiII(1260)':           '-0.780E-3   0.01    None    None    {}'.format(fix_SiII1260_bias_eta),
                       'beta_SiII(1260)':               '0.5        0.      None    None    fixed',
                       'alpha_SiII(1260)':              '1.0        0.      None    None    fixed',
                       'bias_eta_SiIII(1207)':          '-1.716E-3   0.01    None    None    {}'.format(fix_SiIII1207_bias_eta),
                       'beta_SiIII(1207)':              '0.5        0.      None    None    fixed',
                       'alpha_SiIII(1207)':             '1.0        0.      None    None    fixed',
                       'bias_eta_SiII(1193)':           '-1.820E-3   0.01    None    None    {}'.format(fix_SiII1193_bias_eta),
                       'beta_SiII(1193)':               '0.5        0.      None    None    fixed',
                       'alpha_SiII(1193)':              '1.0        0.      None    None    fixed',
                       'bias_eta_SiII(1190)':           '-2.288E-3   0.01    None    None    {}'.format(fix_SiII1190_bias_eta),
                       'beta_SiII(1190)':               '0.5        0.      None    None    fixed',
                       'alpha_SiII(1190)':              '1.0        0.      None    None    fixed',
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

    #create config files for doing fits
    for rmin in args.rmin_values:
        for rmax in args.rmax_values:
            for afix in args.afix_values:
                print(' -> making files for rmin={}, rmax={}, afix={}'.format(rmin,rmax,afix))

                if args.fit_lya_auto | args.fit_all:
                    print(' -> -> making lya_auto files')
                    make_lya_auto_fit_files(analysis_dir.fitsdir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_qso_auto | args.fit_all:
                    print(' -> -> making qso_auto files')
                    make_qso_auto_fit_files(analysis_dir.fitsdir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_dla_auto | args.fit_all:
                    print(' -> -> making dla_auto files')
                    make_dla_auto_fit_files(analysis_dir.fitsdir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_aa_auto | args.fit_all:
                    print(' -> -> making lya_aa_auto files')
                    make_lya_aa_auto_fit_files(analysis_dir.fitsdir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_qso_cross | args.fit_all:
                    print(' -> -> making lya_qso_cross files')
                    make_lya_qso_cross_fit_files(analysis_dir.fitsdir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_dla_cross | args.fit_all:
                    print(' -> -> making lya_dla_cross files')
                    make_lya_dla_cross_fit_files(analysis_dir.fitsdir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_qso_dla_cross | args.fit_all:
                    print(' -> -> making qso_dla_cross files')
                    make_qso_dla_cross_fit_files(analysis_dir.fitsdir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_auto__lya_qso_cross | args.fit_all:
                    print(' -> -> making joint lya_auto__lya_qso_cross files')
                    make_lya_auto__lya_qso_cross_fit_files(analysis_dir.fitsdir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_auto__lya_qso_cross__qso_auto | args.fit_all:
                    print(' -> -> making joint lya_auto__lya_qso_cross__qso_auto files')
                    make_lya_auto__lya_qso_cross__qso_auto_fit_files(analysis_dir.fitsdir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_auto__lya_dla_cross | args.fit_all:
                    print(' -> -> making joint lya_auto__lya_dla_cross files')
                    make_lya_auto__lya_dla_cross_fit_files(analysis_dir.fitsdir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_lya_auto__lya_dla_cross__dla_auto | args.fit_all:
                    print(' -> -> making joint lya_auto__lya_dla_cross__dla_auto files')
                    make_lya_auto__lya_dla_cross__dla_auto_fit_files(analysis_dir.fitsdir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)
                if args.fit_qso_auto__qso_dla_cross__dla_auto | args.fit_all:
                    print(' -> -> making joint qso_auto__qso_dla_cross__dla_auto files')
                    make_qso_auto__qso_dla_cross__dla_auto_fit_files(analysis_dir.fitsdir,exp_filepaths,rmin=rmin,rmax=rmax,afix=afix)



"""
