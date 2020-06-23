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

# Set up options
out_path = args.out_dir + '/master_randoms.fits'
cat_path = args.out_dir + '/master.fits'
nz_filename = 'input_files/Nz_qso_130618_2_colore1_hZs.txt'
max_cat_z = 3.79

# Execute
randoms.generate_rnd(factor=args.rand_factor_QSO,out_path=out_path,
    method=args.rand_method_QSO,cat_path=cat_path,
    footprint=args.footprint,nz_filename=nz_filename,min_cat_z=args.min_cat_z,
    max_cat_z=max_cat_z,overwrite=args.overwrite,N_side=args.nside,
    start_MOCKID_rnd=args.rand_mockid_start,seed=args.seed)
