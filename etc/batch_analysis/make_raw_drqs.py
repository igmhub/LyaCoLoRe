#!/usr/bin/env python

import argparse
import fitsio
import numpy as np
import os

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

parser.add_argument('--qq-ref-zcat',
                    type=str,
                    default=None,
                    required=True,
                    help='Reference zcat file to make sure raw drqs have correct downsampling')

parser.add_argument('--overwrite',
                    action='store_true',
                    default=False,
                    help='Overwrite any existing catalog at the output location')

parser.add_argument('--seed',
                    type=int,
                    default=0,
                    required=False,
                    help='Random seed for downsampling randoms appropriately')

args = parser.parse_args()

def master_to_drq(in_path, out_path, randoms=False, zcat=None, nobj=None):

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

    if nobj is not None:
        gen = np.random.default_rng(seed=args.seed)
        w = gen.permutation(np.arange(len(cat['THING_ID'])))/nobj < 1
        for k in cat.keys():
            cat[k] = cat[k][w]

    # save results
    results = fitsio.FITS(out_path,clobber=args.overwrite)
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
    results = fitsio.FITS(out_path,clobber=args.overwrite)
    cols = list(cat.values())
    names = list(cat)
    results.write(cols, names=names, extname='DLACAT')
    results.close()

    return

## Make raw drq
master = os.path.join(args.in_dir,'master.fits')
out = os.path.join(args.out_dir,'drq_qso.fits')
master_to_drq(master,out,randoms=False,zcat=args.qq_ref_zcat)

## Make qso randoms drq
master_rand = os.path.join(args.in_dir,'master_randoms.fits')

# Work out the multiplier used for randoms
with fitsio.FITS(master) as m:
    nobj = len(m[1][:])
    with fitsio.FITS(master_rand) as mr:
        nobj_rand = len(mr[1][:])
        mult = nobj_rand/nobj
with fitsio.FITS(out) as drq:
    nobj_rand_drq = int(len(drq[1][:])*mult)

out_rand = os.path.join(args.randoms_out_dir,'drq_qso_randoms.fits')
master_to_drq(master_rand,out_rand,randoms=True,zcat=args.qq_ref_zcat,nobj=nobj_rand_drq)

## Make raw dla drq
master_dla = os.path.join(args.in_dir,'master_DLA.fits')
out_dla = os.path.join(args.out_dir,'drq_dla.fits')
master_dla_to_drq(master_dla,out_dla,randoms=False,zcat=args.qq_ref_zcat)

## Make raw dla randoms drq
master_dla_rand = os.path.join(args.in_dir,'master_DLA_randoms.fits')
out_dla_rand = os.path.join(args.randoms_out_dir,'drq_dla_randoms.fits')
master_dla_to_drq(master_dla_rand,out_dla_rand,randoms=True,zcat=args.qq_ref_zcat)
