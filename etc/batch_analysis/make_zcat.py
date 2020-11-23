#!/usr/bin/env python

import argparse
import os

from lyacolore import submit_utils
from subprocess import call

parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )

parser.add_argument('--qq-dir',
                    type=str,
                    default=None,
                    required=True,
                    help='Directory to transmission files')

parser.add_argument('--nside',
                    type=int,
                    default=16,
                    required=False,
                    help='Value of nside used when running picca')

parser.add_argument('--desi-env',
                    type=str,
                    default='20.8',
                    required=False,
                    help='Version of DESI env to use')

args = parser.parse_args()

## Make the text of the script
text = '#!/bin/bash -l\n\n'
text += 'source /global/common/software/desi/desi_environment.sh {}\n'.format(args.desi_env)
spec_dir = os.path.join(args.qq_dir,'spectra-{}'.format(args.nside))
zcat = os.path.join(args.qq_dir,'zcat.fits')
text += 'desi_zcatalog -i {} -o {} --fibermap\n'.format(spec_dir,zcat)

## Write it to file
script = os.path.join(args.qq_dir,'make_zcat.sh')
with open(script,'w') as f:
    f.write(text)
submit_utils.make_file_executable(script)

## Execute it
call(script)
