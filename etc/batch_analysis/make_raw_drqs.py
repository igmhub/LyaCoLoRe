#!/usr/bin/env python

import argparse
import fitsio
import healpy as hp
import numpy as np
import os

from lyacolore import submit_utils
from subprocess import call

parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )

parser.add_argument('-i','--in-dir',
                    type=str,
                    default=None,
                    required=True,
                    help='Directory to transmission files')

parser.add_argument('-o','--out-dir',
                    type=str,
                    default=None,
                    required=True,
                    help='Directory to save drqs in')

parser.add_argument('--randoms-out-dir',
                    type=str,
                    default=None,
                    required=True,
                    help='Directory to save random drqs in')

parser.add_argument('--randoms-zmin',
                    type=float,
                    default=1.7,
                    required=False,
                    help='Min value of z for QSOs (applied to randoms)')

parser.add_argument('--randoms-downsampling',
                    type=float,
                    default=0.4,
                    required=False,
                    help='Proportion by which to downsample (applied to randoms)')

parser.add_argument('--qq-ref-zcat',
                    type=str,
                    default=None,
                    required=True,
                    help='Reference zcat file to make sure raw drqs have correct downsampling')

parser.add_argument('--seed',
                    type=int,
                    default=0,
                    required=False,
                    help='Random seed for downsampling randoms appropriately')

parser.add_argument('--desi-env',
                    type=str,
                    default='master',
                    required=False,
                    help='DESI env to use for making randoms have right footprint')

args = parser.parse_args()

def master_to_drq(in_path, out_path, randoms=False, zcat=None, randoms_downsampling=None, randoms_zmin=None):

    from_desi_key_to_picca_key = {
        'RA': 'RA',
        'DEC': 'DEC',
        'Z': 'Z_QSO_RSD',
        'THING_ID': 'MOCKID',
        'PLATE': 'MOCKID',
        'MJD': 'MOCKID',
        'FIBERID': 'MOCKID'
    }
    if randoms:
        from_desi_key_to_picca_key['Z'] = 'Z'
    # read catalogue
    cat = {}
    hdul = fitsio.FITS(in_path)
    for key, value in from_desi_key_to_picca_key.items():
        cat[key] = hdul['CATALOG'][value][:]
    hdul.close()
    print(("INFO: Found {} quasars").format( np.unique(cat['THING_ID']).size))

    # sort by THING_ID
    w = np.argsort(cat['THING_ID'])
    for key in cat:
        cat[key] = cat[key][w]

    for key in ['RA', 'DEC']:
        cat[key] = cat[key].astype('float64')

    if (zcat is not None) and (not randoms):
        with fitsio.FITS(zcat) as zc:
            thingid_zcat = zc[1][:]['TARGETID']
        w = np.in1d(cat['THING_ID'],thingid_zcat)
        print('INFO: zcat contains {} quasars, of which {} found in master.'.format(len(thingid_zcat),w.sum()))
        print('INFO: Reducing to this set now.')
        for k in cat.keys():
            cat[k] = cat[k][w]

    elif randoms:
        w = np.ones(cat['THING_ID'].shape).astype(bool)
        if randoms_zmin is not None:
            w &= (cat['Z']>=randoms_zmin)
        if randoms_downsampling is not None:
            gen = np.random.default_rng(seed=args.seed)
            w &= (gen.permutation(np.arange(len(cat['THING_ID'])))/len(cat['THING_ID']) < randoms_downsampling)

    # save results
    results = fitsio.FITS(out_path, 'rw', clobber=True)
    cols = list(cat.values())
    names = list(cat)
    results.write(cols, names=names, extname='CAT')
    results.close()

    return

def master_dla_to_drq(in_path, out_path, randoms=False, zcat=None):

    from_desi_key_to_picca_key = {
        'RA': 'RA',
        'DEC': 'DEC',
        'Z': 'Z_DLA_RSD',
        'ZQSO': 'Z_QSO_RSD',
        'NHI': 'N_HI_DLA',
        'THING_ID': 'MOCKID',
        'DLAID': 'DLAID',
        'PLATE': 'MOCKID',
        'MJD': 'MOCKID',
        'FIBERID': 'MOCKID',
    }
    if randoms:
        from_desi_key_to_picca_key['Z'] = 'Z_DLA'
        from_desi_key_to_picca_key.pop('NHI')
    # read catalogue
    cat = {}
    hdul = fitsio.FITS(in_path)
    for key, value in from_desi_key_to_picca_key.items():
        cat[key] = hdul['DLACAT'][value][:]
    hdul.close()
    print(("INFO: Found {} DLA from {} "
               "quasars").format(cat['Z'].size,
                                 np.unique(cat['THING_ID']).size))
    # sort by THING_ID
    w = np.argsort(cat['THING_ID'])
    for key in cat:
        cat[key] = cat[key][w]

    for key in ['RA', 'DEC']:
        cat[key] = cat[key].astype('float64')

    if zcat is not None:
        with fitsio.FITS(zcat) as zc:
            thingid_zcat = zc[1][:]['TARGETID']
        w = np.in1d(cat['THING_ID'],thingid_zcat)
        print('INFO: zcat contains {} quasars, of which {} have DLAs.'.format(len(thingid_zcat),w.sum()))
        print('INFO: Reducing to this set now.')
        for k in cat.keys():
            cat[k] = cat[k][w]

    # save results
    results = fitsio.FITS(out_path, 'rw', clobber=True)
    cols = list(cat.values())
    names = list(cat)
    results.write(cols, names=names, extname='DLACAT')
    results.close()

    return

## Make raw drq
print('INFO: Making raw drq')
master = os.path.join(args.in_dir,'master.fits')
out = os.path.join(args.out_dir,'drq_qso.fits')
master_to_drq(master,out,randoms=False,zcat=args.qq_ref_zcat)
print('')

## Make qso randoms drq
print('INFO: Making randoms drq')
master_rand = os.path.join(args.in_dir,'master_randoms.fits')
out_rand = os.path.join(args.randoms_out_dir,'drq_qso_randoms.fits')
master_to_drq(master_rand,out_rand,randoms=True,zcat=args.qq_ref_zcat,randoms_downsampling=args.randoms_downsampling,randoms_zmin=args.randoms_zmin)
print('')

## Make raw dla drq
print('INFO: Making raw dla drq')
master_dla = os.path.join(args.in_dir,'master_DLA.fits')
out_dla = os.path.join(args.out_dir,'drq_dla.fits')
master_dla_to_drq(master_dla,out_dla,randoms=False,zcat=args.qq_ref_zcat)
print('')

## Make raw dla randoms drq
print('INFO: Making dla randoms drq')
master_dla_rand = os.path.join(args.in_dir,'master_DLA_randoms.fits')
out_dla_rand = os.path.join(args.randoms_out_dir,'drq_dla_randoms.fits')
master_dla_to_drq(master_dla_rand,out_dla_rand,randoms=True,zcat=args.qq_ref_zcat)
print('')

## Make script to fix randoms drq footprint.
print('INFO: Correcting footprint of randoms drq')

# Make the text of the script
text = '#!/bin/bash -l\n\n'
text += 'source /global/common/software/desi/desi_environment.sh {}\n'.format(args.desi_env)
text += '/global/homes/j/jfarr/Projects/LyaCoLoRe/etc/batch_analysis/apply_randoms_footprint.py --randoms-dir {}'.format(args.randoms_out_dir)

# Write it to file
script = os.path.join(args.randoms_out_dir,'apply_randoms_footprint.sh')
with open(script,'w') as f:
    f.write(text)
submit_utils.make_file_executable(script)

# Execute it
call(script)
