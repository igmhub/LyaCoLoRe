#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import argparse
from scipy.interpolate import interp1d
import h5py

ver = 'stack'
fit_types = ['lya_auto','lya_qso_cross','lya_auto__lya_qso_cross','lya_dla_cross','lya_auto__lya_dla_cross','qso_auto']
rmins = [20.,40.0]
pars = ['ap', 'at', 'bias_eta_LYA', 'beta_LYA', 'beta_QSO', 'beta_DLA']

nchar = np.max([len(par) for par in pars])
print(' ')
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

        f = ff['best fit'].attrs['growth_rate'][0]
        try:
            b_LYA = ff['best fit'].attrs['bias_eta_LYA'][0]*f/ff['best fit'].attrs['beta_LYA'][0]
            b_LYA_err = b_LYA * np.sqrt((ff['best fit'].attrs['bias_eta_LYA'][1]/ff['best fit'].attrs['bias_eta_LYA'][0])**2 + (ff['best fit'].attrs['beta_LYA'][1]/ff['best fit'].attrs['beta_LYA'][0])**2)
            print('{:12s} = {:2.3f} \pm {:1.4f}'.format('b_F',b_LYA,b_LYA_err))
        except:
            pass
        try:    
            b_QSO = f/ff['best fit'].attrs['beta_QSO'][0]
            b_QSO_err = b_QSO * (ff['best fit'].attrs['beta_QSO'][1]/ff['best fit'].attrs['beta_QSO'][0])
            print('{:12s} = {:2.3f} \pm {:1.4f}'.format('b_QSO',b_QSO,b_QSO_err))
        except:
            pass
        try:    
            b_DLA = f/ff['best fit'].attrs['beta_DLA'][0]
            b_DLA_err = b_DLA * (ff['best fit'].attrs['beta_DLA'][1]/ff['best fit'].attrs['beta_DLA'][0])
            print('{:12s} = {:2.3f} \pm {:1.4f}'.format('b_DLA',b_DLA,b_DLA_err))
        except:
            pass
        print(' ')
