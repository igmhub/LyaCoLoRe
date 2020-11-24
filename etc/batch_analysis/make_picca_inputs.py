#!/usr/bin/env python

import argparse
import os

from lyacolore import submit_utils
from subprocess import call

parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )

parser.add_argument('--raw-dir',
                    type=str,
                    default=None,
                    required=True,
                    help='Raw directory for a given realisation')

parser.add_argument('--qq-basedir',
                    type=str,
                    default=None,
                    required=True,
                    help='Directory that contains all of the qq runs for the given realisation')

parser.add_argument('--qq-runs',
                    type=str,
                    nargs='*',
                    default=None,
                    required=True,
                    help='Directory names of the qq runs')

parser.add_argument('--picca-basedir',
                    type=str,
                    default=None,
                    required=True,
                    help='Directory that contains all of the picca outputs for the given realisation')

args = parser.parse_args()

## For each of a sequence of qq runs:
for qq_run in args.qq_runs:

    ## Point to the qq output.
    qq_dir = os.path.join(args.qq_basedir,qq_run)

    ## Make the directory structure for our outputs.
    analysis_dir = submit_utils.AnalysisDir(args.picca_basedir,'desi-raw')

    ## Make the zcat
    submit_utils.make_permission_group_desi(qq_dir)
    command = '/global/homes/j/jfarr/Projects/LyaCoLoRe/etc/batch_analysis/make_zcat.py --qq-dir {}'.format(qq_dir)
    call(command.split(' '))

    ## Make the drq
    submit_utils.make_permission_group_desi(qq_dir)
    command = '/global/homes/j/jfarr/Projects/LyaCoLoRe/etc/batch_analysis/make_drqs.py --in-dir {} --out-dir {}'.format(qq_dir,analysis_dir.datadir)
    call(command.split(' '))

    ## Submit job to run the deltas
    drq = os.path.join(qq_out_data_dir,'drq_qso.fits')
    in_dir = os.path.join(qq_out_data_dir,'spectra-16')
    command = '/global/homes/j/jfarr/Projects/LyaCoLoRe/etc/batch_analysis/run_make_deltas.py --out-dir {} --drq {} --in-dir {}'.format(analysis_dir.datadir,drq,in_dir)
    call(command.split(' '))


## Point to the raw output.
raw_in_dir = args.raw_dir

## Make the directory structure for our outputs.
analysis_dir = submit_utils.AnalysisDir(args.picca_basedir,'desi-raw')

## Point to the picca directory for the randoms output, make sure it exists.
randoms_dir = os.path.join(args.picca_basedir,'randoms')
if not os.path.isdir(randoms_dir):
    os.mkdir(randoms_dir)
submit_utils.make_permission_group_desi(randoms_dir)

## Make the raw drqs
qq_ref_zcat = os.path.join(args.qq_basedir,args.qq_runs[0],'zcat.fits')
command = '/global/homes/j/jfarr/Projects/LyaCoLoRe/etc/batch_analysis/make_raw_drqs.py --in-dir {} --out-dir {} --randoms-out-dir {} --qq-ref-zcat {}'.format(raw_in_dir,analysis_dir.datadir,randoms_dir,qq_ref_zcat)
call(command.split(' '))

## Submit job to run the raw deltas
drq = os.path.join(raw_out_data_dir,'drq_qso.fits')
in_dir = args.raw_dir
command = '/global/homes/j/jfarr/Projects/LyaCoLoRe/etc/batch_analysis/run_make_raw_deltas.py --out-dir {} --drq {} --in-dir {}'.format(analysis_dir.datadir,drq,in_dir)
call(command.split(' '))
