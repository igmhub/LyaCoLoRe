#!/usr/bin/env python

import sys

from lyacolore import parse, randoms

################################################################################

"""
Set up the file locations and filename structures.
Also define option preferences.
"""

args = parse.get_args(sys.argv)

################################################################################

#Set up options
out_path = args.out_dir + '/master_DLA_randoms.fits'
DLA_cat_path = args.out_dir + '/master_DLA.fits'
QSO_cat_path = args.out_dir + '/master.fits'

# Execute
randoms.generate_rnd_dla(factor=args.rand_factor_DLA,out_path=out_path,
    DLA_cat_path=DLA_cat_path,QSO_cat_path=QSO_cat_path,
    footprint=args.footprint,lambda_min=args.transmission_lambda_min,
    lambda_max=args.transmission_lambda_max,NHI_min=args.DLA_NHI_min,
    NHI_max=args.DLA_NHI_max,overwrite=args.overwrite,N_side=args.nside,
    add_NHI=args.add_rand_DLA_NHI,method=args.rand_method_DLA,
    start_DLAID_rnd=args.rand_dlaid_start,seed=args.seed)
