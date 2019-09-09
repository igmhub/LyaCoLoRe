#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import argparse
from scipy.interpolate import interp1d
import h5py

location = '/project/projectdirs/desi/users/jfarr/LyaCoLoRe_paper/analysis/correlation_functions/stack/fits/lya_aa_auto/'
result_filepath = glob.glob(location+'*40*.h5')
ff = h5py.File(result_filepath,'r')
fit = {}

pars = ['ap', 'at', 'bias_eta_LYA', 'beta_LYA', 'beta_QSO']
for par in pars:
    print('{} = {:1.3f} \pm {:1.4f}'.format(par,ff['best fit'].attrs[par][0],ff['best fit'].attrs[par][1]))

f = ff['best fit'].attrs['growth_rate']
b_LYA = ff['best fit'].attrs['bias_eta_LYA']*f/ff['best fit'].attrs['beta_LYA']
print('{} = {:1.3f} \pm {:1.4f}'.format('b_F',ff['best fit'].attrs[par][0],ff['best fit'].attrs[par][1]))
b_QSO = f/ff['best fit'].attrs['beta_QSO']
print('{} = {:1.3f} \pm {:1.4f}'.format('b_QSO',ff['best fit'].attrs[par][0],ff['best fit'].attrs[par][1]))
