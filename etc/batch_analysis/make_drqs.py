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
                    help='Directory containing zcat')

parser.add_argument('-o','--out-dir',
                    type=str,
                    default=None,
                    required=True,
                    help='Directory to save drqs to')

parser.add_argument('--overwrite',
                    action='store_true',
                    default=False,
                    help='Overwrite any existing catalog at the output location')

args = parser.parse_args()

def zcat_to_drq(in_path, out_path):

    from_desi_key_to_picca_key = {
        'RA': 'RA',
        'DEC': 'DEC',
        'Z': 'Z',
        'THING_ID': 'TARGETID',
        'PLATE': 'TARGETID',
        'MJD': 'TARGETID',
        'FIBERID': 'TARGETID'
    }

    # read catalogue
    cat = {}
    hdul = fitsio.FITS(in_path)
    for key, value in from_desi_key_to_picca_key.items():
        cat[key] = hdul['ZCATALOG'][value][:]
    hdul.close()
    print(("INFO: Found {} quasars").format( np.unique(cat['THING_ID']).size))

    # sort by THING_ID
    w = np.argsort(cat['THING_ID'])
    for key in cat:
        cat[key] = cat[key][w]

    for key in ['RA', 'DEC']:
        cat[key] = cat[key].astype('float64')

    # save results
    results = fitsio.FITS(out_path, 'rw', clobber=args.overwrite)
    cols = list(cat.values())
    names = list(cat)
    results.write(cols, names=names, extname='CAT')
    results.close()

## Make drq
zcat = os.path.join(args.in_dir,'zcat.fits')
drq = os.path.join(args.out_dir,'drq_qso.fits')
zcat_to_drq(zcat,drq)
