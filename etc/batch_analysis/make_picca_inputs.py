#!/usr/bin/env python

import argparse
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

parser.add_argument('--picca-basedir',
                    type=str,
                    default=None,
                    required=True,
                    help='Directory that contains all of the picca outputs for the given realisation')

args = parser.parse_args()

## Make the 1.0 zcat
## This needs to be done before the raw drqs and deltas
qq_dir = os.path.join(args.qq_basedir,'desi-1.0-4')
if not os.path.isdir(qq_dir):
    os.mkdir(qq_dir)
command = 'make_zcats.py --qq-dir {}'.format(qq_dir)
call(command)



## Make the raw drqs
raw_in_dir = args.raw_dir
raw_out_dir = os.path.join(args.picca_basedir,'desi-raw')
if not os.path.isdir(raw_out_dir):
    os.mkdir(raw_out_dir)
randoms_dir = os.path.join(args.picca_basedir,'randoms')
if not os.path.isdir(randoms_dir):
    os.mkdir(randoms_dir)
qq_ref_zcat = os.path.join(qq_dir,'zcat.fits')
command = 'make_raw_drqs.py --in-dir {} --out-dir {} --randoms-out-dir {} --qq-ref-zcat {}'.format(raw_in_dir,raw_out_dir,randoms_dir,qq_ref_zcat)
call(command)

## Submit job to run the raw deltas
python_script = os.path.join(args.picca_basedir,'desi-raw/run_picca_deltas.py')
slurm_script = os.path.join(args.picca_basedir,'desi-raw/run_picca_deltas.sh')
slurm_hours = 2.
out_dir = os.path.join(args.picca_basedir,'desi-raw/deltas')
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
drq = os.path.join(args.raw_dir,'drq.fits')
in_dir = args.raw_dir
command = 'run_make_raw_deltas.py --python-script {} --slurm-script {} --slurm-hours {} --out-dir {} --drq {} --in-dir {}'.format(python_script,slurm_script,slurm_hours,out_dir,drq,in_dir)



## Make the 1.0 drq
qq_dir = os.path.join(args.qq_basedir,'desi-1.0-4')
if not os.path.isdir(qq_dir):
    os.mkdir(qq_dir)
command = 'make_drqs.py --in-dir {}'.format(qq_dir)
call(command)

## Submit job to run the 1.0 deltas
slurm_script = os.path.join(args.picca_basedir,'desi-1.0-4/run_picca_deltas.sh')
slurm_hours = 2.
out_dir = os.path.join(args.picca_basedir,'desi-1.0-4/deltas')
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
drq = os.path.join(args.qq_basedir,'desi-1.0-4/drq.fits')
in_dir = os.path.join(args.qq_basedir,'desi-1.0-4')
command = 'run_make_deltas.py --slurm-script {} --slurm-hours {} --out-dir {} --drq {} --in-dir {}'.format(slurm_script,slurm_hours,out_dir,drq,in_dir)
call(command)
