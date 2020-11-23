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

    ## Make the zcat
    qq_dir = os.path.join(args.qq_basedir,qq_run)
    if not os.path.isdir(qq_dir):
        os.mkdir(qq_dir)
    submit_utils.make_permission_group_desi(qq_dir)
    command = 'make_zcats.py --qq-dir {}'.format(qq_dir)
    call(command.split(' '))



## Make the raw drqs
raw_in_dir = args.raw_dir
raw_out_dir = os.path.join(args.picca_basedir,'desi-raw')
if not os.path.isdir(raw_out_dir):
    os.mkdir(raw_out_dir)
submit_utils.make_permission_group_desi(qq_dir)
randoms_dir = os.path.join(args.picca_basedir,'randoms')
if not os.path.isdir(randoms_dir):
    os.mkdir(randoms_dir)
submit_utils.make_permission_group_desi(randoms_dir)
qq_ref_zcat = os.path.join(qq_dir,'zcat.fits')
command = 'make_raw_drqs.py --in-dir {} --out-dir {} --randoms-out-dir {} --qq-ref-zcat {}'.format(raw_in_dir,raw_out_dir,randoms_dir,qq_ref_zcat)
call(command)

## Submit job to run the raw deltas
python_script = os.path.join(args.picca_basedir,'desi-raw/run_picca_deltas.py')
slurm_script = os.path.join(args.picca_basedir,'desi-raw/run_picca_deltas.sl')
slurm_hours = 0.5
slurm_queue = 'debug'
out_dir = os.path.join(args.picca_basedir,'desi-raw/deltas')
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
submit_utils.make_permission_group_desi(out_dir)
drq = os.path.join(args.picca_basedir,'drq_qso.fits')
in_dir = args.raw_dir
command = 'run_make_raw_deltas.py --python-script {} --slurm-script {} --slurm-hours {} --slurm-queue {} --out-dir {} --drq {} --in-dir {}'.format(python_script,slurm_script,slurm_hours,slurm_queue,out_dir,drq,in_dir)



## For each of a sequence of qq runs:
for qq_run in args.qq_runs:

    ## Make the drq
    qq_dir = os.path.join(args.qq_basedir,qq_run)
    qq_out_dir = os.path.join(args.picca_basedir,qq_run)
    if not os.path.isdir(qq_dir):
        os.mkdir(qq_dir)
    submit_utils.make_permission_group_desi(qq_dir)
    command = 'make_drqs.py --in-dir {} --out-dir {}'.format(qq_dir,qq_out_dir)
    call(command)

    ## Submit job to run the deltas
    slurm_script = os.path.join(args.picca_basedir,qq_run,'run_picca_deltas.sl')
    slurm_hours = 2.
    slurm_queue = 'regular'
    out_dir = os.path.join(args.picca_basedir,qq_run)
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    deltas_dir = os.path.join(args.picca_basedir,qq_run,'deltas/')
    if not os.path.isdir(deltas_dir):
        os.mkdir(deltas_dir)
    submit_utils.make_permission_group_desi(out_dir)
    drq = os.path.join(args.picca_basedir,qq_run,'drq_qso.fits')
    in_dir = os.path.join(args.qq_basedir,qq_run)
    command = 'run_make_deltas.py --slurm-script {} --slurm-hours {} --slurm-queue {} --out-dir {} --drq {} --in-dir {}'.format(slurm_script,slurm_hours,slurm_queue,out_dir,drq,in_dir)
    call(command)
