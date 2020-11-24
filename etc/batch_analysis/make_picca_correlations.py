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

fid_Om = 0.314569514863487
fid_Or = 7.97505418919554e-05

randoms_dir = os.path.join(args.picca_basedir,'randoms')

## Make the directory structure for our outputs.
raw_analysis_dir = submit_utils.AnalysisDir(args.picca_basedir,'desi-raw')

## For each of a sequence of qq runs:
for qq_run in args.qq_runs:

    ## Point to the qq output.
    qq_dir = os.path.join(args.qq_basedir,qq_run)

    ## Make the directory structure for our outputs.
    analysis_dir = submit_utils.AnalysisDir(args.picca_basedir,qq_run)

    ## Setting up the correlation jobs
    print('INFO: Setting up correlations from picca data in {}'.format(analysis_dir.datadir))
    command = '/global/homes/j/jfarr/Projects/LyaCoLoRe/etc/batch_analysis/run_picca_correlations.py'
    command += ' --deltas-dir {}'.format(analysis_dir.deltadir)
    command += ' --drq-qso {}'.format(os.path.join(analysis_dir.datadir,'drq_qso.fits'))
    command += ' --drq-dla {}'.format(os.path.join(raw_analysis_dir.datadir,'drq_dla.fits'))
    command += ' --drq-qso-randoms {}'.format(os.path.join(randoms_dir,'drq_qso_randoms.fits'))
    command += ' --drq-dla-randoms {}'.format(os.path.join(randoms_dir,'drq_dla_randoms.fits'))
    command += ' --corr-dir {}'.format(analysis_dirr.corrdir)
    command += ' --fid-Om {}'.format(fid_Om)
    command += ' --fid-Or {}'.format(fid_Or)
    command += ' --run-lya-auto'
    command += ' --run-lya-qso-cross'
    command += ' --no-submit'
    command += '\n'
    call(command.split(' '))
    print('')


## Setting up the correlation jobs
print('INFO: Setting up correlations from picca data in {}'.format(raw_analysis_dir.datadir))
command = '/global/homes/j/jfarr/Projects/LyaCoLoRe/etc/batch_analysis/run_picca_correlations.py'
command += ' --deltas-dir {}'.format(raw_analysis_dir.deltadir)
command += ' --drq-qso {}'.format(os.path.join(raw_analysis_dir.datadir,'drq_qso.fits'))
command += ' --drq-dla {}'.format(os.path.join(raw_analysis_dir.datadir,'drq_dla.fits'))
command += ' --drq-qso-randoms {}'.format(os.path.join(randoms_dir,'drq_qso_randoms.fits'))
command += ' --drq-dla-randoms {}'.format(os.path.join(randoms_dir,'drq_dla_randoms.fits'))
command += ' --corr-dir {}'.format(raw_analysis_dir.corrdir)
command += ' --fid-Om {}'.format(fid_Om)
command += ' --fid-Or {}'.format(fid_Or)
command += ' --run-lya-auto'
command += ' --run-lya-qso-cross'
command += ' --no-project'
command += ' --no-submit'
command += '\n'
call(command.split(' '))
print('')
