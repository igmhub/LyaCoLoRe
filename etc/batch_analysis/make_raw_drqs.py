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

args = parser.parse_args()

def master_to_drq(in_path, out_path, randoms=False, zcat=None):

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

    if zcat is not None:
        with fits.open(zcat) as zc:
            thingid_zcat = zc[1].data['TARGETID']
        w = np.in1d(cat['THING_ID'],thingid_zcat)
        print('INFO: zcat contains {} quasars, of which {} found in master.'.format(len(thingid_zcat),w.sum()))
        print('INFO: Reducing to this set now.')
        for k in cat.keys():
            cat[k] = cat[k][w]

    # save results
    results = fitsio.FITS(out_path, 'rw', clobber=args.overwrite)
    cols = list(cat.values())
    names = list(cat)
    results.write(cols, names=names, extname='CAT')
    results.close()

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
    userprint(("INFO: Found {} DLA from {} "
               "quasars").format(cat['Z'].size,
                                 np.unique(cat['THING_ID']).size))
    # sort by THING_ID
    w = np.argsort(cat['THING_ID'])
    for key in cat:
        cat[key] = cat[key][w]

    for key in ['RA', 'DEC']:
        cat[key] = cat[key].astype('float64')

    if zcat is not None:
        with fits.open(zcat) as zc:
            thingid_zcat = zc[1].data['TARGETID']
        w = np.in1d(cat['THING_ID'],thingid_zcat)
        print('INFO: zcat contains {} quasars, of which {} found in master.'.format(len(thingid_zcat),w.sum()))
        print('INFO: Reducing to this set now.')
        for k in cat.keys():
            cat[k] = cat[k][w]

    # save results
    results = fitsio.FITS(out_path, 'rw', clobber=args.overwrite)
    cols = list(cat.values())
    names = list(cat)
    results.write(cols, names=names, extname='DLACAT')
    results.close()

## Make raw drq
master = os.path.join(args.in_dir,'master.fits')
out = os.path.join(args.out_dir,'drq.fits')
master_to_drq(master,out,randoms=False,zcat=args.qq_ref_zcat)

## Make qso randoms drq
master = os.path.join(args.in_dir,'master_randoms.fits')
out = os.path.join(args.randoms_out_dir,'drq_randoms.fits')
master_to_drq(master,out,randoms=True,zcat=args.qq_ref_zcat)

## Make raw dla drq
master = os.path.join(args.in_dir,'master_dla.fits')
out = os.path.join(args.out_dir,'drq_dla.fits')
master_dla_to_drq(master,out,randoms=False,zcat=args.qq_ref_zcat)

## Make raw dla randoms drq
master = os.path.join(args.in_dir,'master_dla_randoms.fits')
out = os.path.join(args.randoms_out_dir,'drq_dla_randoms.fits')
master_dla_to_drq(master,out,randoms=True,zcat=args.qq_ref_zcat)
