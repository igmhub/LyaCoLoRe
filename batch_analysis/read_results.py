#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import argparse
from scipy.interpolate import interp1d
import h5py

ver = 'stack'
fit_types = ['lya_auto','lya_qso_cross','lya_auto__lya_qso_cross','lya_dla_cross','lya_auto__lya_dla_cross','qso_auto','lya_aa_auto']
rmins = [20.,40.0]
pars = ['ap', 'at']
tracers = ['LYA','QSO','DLA','SiII(1190)','SiII(1193)','SiIII(1207)','SiII(1260)']

nchar = np.max([len(par) for par in pars])
print(' ')
print('-'*80)
for fit_type in fit_types:
    print(fit_type)
    for rmin in rmins:
        print('rmin =',rmin)
        location = '/project/projectdirs/desi/users/jfarr/LyaCoLoRe_paper/analysis/correlation_functions/{}/fits/{}/'.format(ver,fit_type)
        result_filepath = location+'result_{}_rmin{}_rmax160.0_afree.h5'.format(fit_type,rmin)
        ff = h5py.File(result_filepath,'r')
        fit = {}

        for par in pars:
            try:
                print('{:12s} = {:2.3f} \pm {:1.4f}'.format(par,ff['best fit'].attrs[par][0],ff['best fit'].attrs[par][1]))
            except:
                print('{:12s} not found'.format(par))

        for t in tracers:
            try:
                print(t)
                b_eta = ff['best fit'].attrs['bias_eta_{}'.format(t)][0]
                b_eta_err = ff['best fit'].attrs['bias_eta_{}'.format(t)][1]
                print('{:12s} = {:2.3e} \pm {:1.4e}'.format('bias_eta',b_eta,b_eta_err))

                beta = ff['best fit'].attrs['beta_{}'.format(t)][0]
                beta_err = ff['best fit'].attrs['beta_{}'.format(t)][1]
                print('{:12s} = {:2.3e} \pm {:1.4e}'.format('beta',beta,beta_err))

                b = b_eta*f/beta
                b_err = b * np.sqrt((b_eta_err/b_eta)**2 + (beta_err/beta)**2)
                print('{:12s} = {:2.3e} \pm {:1.4e}'.format('bias',b,b_err))
                print(' ')
            except:
                pass

        print(' ')
        print('-'*80)
