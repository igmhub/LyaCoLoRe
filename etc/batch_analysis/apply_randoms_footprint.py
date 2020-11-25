#!/usr/bin/env python

import argparse
import fitsio
import os

from desimodel.footprint import load_tiles, is_point_in_desi

parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )

parser.add_argument('--randoms-dir',
                    type=str,
                    default=None,
                    required=True,
                    help='Directory containing randoms')

parser.add_argument('--randoms-name',
                    type=str,
                    default='drq_qso_randoms.fits',
                    required=False,
                    help='QSO randoms filename')

args = parser.parse_args()

randoms_path = os.path.join(args.randoms_dir,args.randoms_name)
with fitsio.FITS(randoms_path,'rw') as h:
    tiles = load_tiles()
    w = is_point_in_desi(tiles,h[1][:]['RA'],h[1][:]['DEC'])
    cat = h[1][w]

h = fitsio.FITS(randoms_path,'rw',clobber=True)
h.write(cat, extname='CAT')
