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
                    required=False,
                    help='Raw directory for a given realisation')

parser.add_argument('--qq-basedir',
                    type=str,
                    default=None,
                    required=False,
                    help='Directory that contains all of the qq runs for the given realisation')

parser.add_argument('--qq-runs',
                    type=str,
                    nargs='*',
                    default=None,
                    required=False,
                    help='Directory names of the qq runs')

parser.add_argument('--picca-basedir',
                    type=str,
                    default=None,
                    required=True,
                    help='Directory that contains all of the picca outputs for the given realisation')

args = parser.parse_args()

## For each of a sequence of qq runs:
if (args.qq_basedir is not None) and (args.qq_runs is not None):
    for qq_run in args.qq_runs:

        ## Point to the qq output.
        qq_dir = os.path.join(args.qq_basedir,qq_run)

        ## Make the directory structure for our outputs.
        analysis_dir = submit_utils.AnalysisDir(args.picca_basedir,qq_run)

        ## Make the zcat
        print('INFO: Making zcat for quickquasars output in {}'.format(qq_dir))
        command = '/global/homes/j/jfarr/Projects/LyaCoLoRe/etc/batch_analysis/make_zcat.py --qq-dir {}'.format(qq_dir)
        call(command.split(' '))
        print('')

        ## Make the drq
        print('INFO: Making drqs for quickquasars output in {}'.format(qq_dir))
        command = '/global/homes/j/jfarr/Projects/LyaCoLoRe/etc/batch_analysis/make_drqs.py --in-dir {} --out-dir {}'.format(qq_dir,analysis_dir.datadir)
        call(command.split(' '))
        print('')

        ## Submit job to run the deltas
        print('INFO: Submitting job to make deltas for quickquasars output in {}'.format(qq_dir))
        drq = os.path.join(analysis_dir.datadir,'drq_qso.fits')
        in_dir = os.path.join(qq_dir,'spectra-16')
        command = '/global/homes/j/jfarr/Projects/LyaCoLoRe/etc/batch_analysis/run_make_deltas.py --out-dir {} --drq {} --in-dir {}'.format(analysis_dir.datadir,drq,in_dir)
        call(command.split(' '))
        print('')



if args.raw_dir is not None:
    ## Make the directory structure for our outputs.
    analysis_dir = submit_utils.AnalysisDir(args.picca_basedir,'desi-raw')

    ## Point to the picca directory for the randoms output, make sure it exists.
    randoms_dir = os.path.join(args.picca_basedir,'randoms')
    if not os.path.isdir(randoms_dir):
        os.mkdir(randoms_dir)
    submit_utils.make_permission_group_desi(randoms_dir)

    ## Make the raw drqs
    print('INFO: Making drqs for raw output in {}'.format(args.raw_dir))
    qq_ref_zcat = os.path.join(args.qq_basedir,args.qq_runs[0],'zcat.fits')
    command = '/global/homes/j/jfarr/Projects/LyaCoLoRe/etc/batch_analysis/make_raw_drqs.py --in-dir {} --out-dir {} --randoms-out-dir {} --qq-ref-zcat {}'.format(args.raw_dir,analysis_dir.datadir,randoms_dir,qq_ref_zcat)
    call(command.split(' '))
    print('')

    ## Submit job to run the raw deltas
    print('INFO: Submitting job to make deltas for raw output in {}'.format(args.raw_dir))
    drq = os.path.join(analysis_dir.datadir,'drq_qso.fits')
    in_dir = args.raw_dir
    command = '/global/homes/j/jfarr/Projects/LyaCoLoRe/etc/batch_analysis/run_make_raw_deltas.py --out-dir {} --drq {} --in-dir {}'.format(analysis_dir.datadir,drq,in_dir)
    call(command.split(' '))
    print('')
