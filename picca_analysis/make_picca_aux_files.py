#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#locations and names

parser.add_argument('--base-dir', type = str, default = None, required=True,
                    help = 'base directory')

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

#Fit variable

parser.add_argument('--rmin-values', type = float, default = 0., required=True,
                    help = 'minimum separation', nargs='*')

parser.add_argument('--rmax-values', type = float, default = 160., required=True,
                    help = 'maximum separation', nargs='*')

parser.add_argument('--afix-values', type = str, default = None, required=True,
                    help = 'maximum separation', nargs='*')


args = parser.parse_args()


#Calculate the effective redshift.
h = fits.open(args.base_dir+'/'+args.corr_exp_filename)
R = np.sqrt(h[1].data['RP']**2 + h[1].data['RT']**2)
cells = (R > 80.) * (R < 120.)
zeff = np.average(h[1].data['Z'][cells],weights=h[1].data['NB'][cells])
h.close()

#Calculate f. Assume Or=0 for now.
# TODO: include the input Om and Or from picca.
Om_z0=0.3147
Om = Om_z0 * ((1+zeff)**3) / (Om_z0 * ((1+zeff)**3) + 1 - Om_z0)
Ol = 1 - Om
f = (Om**0.6) + (Ol/70.)*(1 + Om/2.)

#Calculate beta_QSO by interpolating the input bias.
b_of_z = np.loadtxt('/global/homes/j/jfarr/Projects/LyaCoLoRe/input_files/Bz_qso_G18.txt')
b_QSO = np.interp(zeff,b_of_z[:,0],b_of_z[:,1])
b_QSO = 2.
beta_QSO = f/b_QSO

#Write parameter file.
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
file = open(args.base_dir+args.parameter_filename,'w')
file.write(parameter_file_text)
file.close()

#create config files for doing fits
for rmin in args.rmin_values:
    for rmax in args.rmax_values:
        for afix in args.afix_values:

            suffix = '{}r_a{}'.format(int(rmin),afix)
            config_filename = args.base_dir+'config_{}_{}.ini'.format(args.corr_type,suffix)
            chi2_filename = args.base_dir+'chi2_{}.ini'.format(suffix)

            chi2_text = ''
            chi2_text += '#!/bin/bash -l\n\n'
            chi2_text += '[data sets]\n'
            chi2_text += 'zeff = {}\n'.format(zeff)
            chi2_text += 'ini files = {}\n\n'.format(config_filename)
            chi2_text += '[cosmo-fit type]\n'
            chi2_text += 'cosmo fit func = ap_at\n\n'
            chi2_text += '[output]\n'
            chi2_text += 'filename = ./result_{}.h5\n\n'.format(suffix)
            chi2_text += '[fiducial]\n'
            chi2_text += 'filename = /PlanckDR12/PlanckDR12.fits\n'
            file = open(chi2_filename,'w')
            file.write(chi2_text)
            file.close()

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
            config_text += 'per_sigma_smooth = 3.78 0.4 None None fixed\n'
            file = open(config_filename,'w')
            file.write(config_text)
            file.close()
