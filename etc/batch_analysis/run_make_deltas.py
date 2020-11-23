#!/usr/bin/env python

import argparse
import os

from subprocess import call
from lyacolore import submit_utils

parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )

parser.add_argument('--slurm-script',
                    type=str,
                    default=None,
                    required=True,
                    help='Output slurm script')

parser.add_argument('--slurm-hours',
                    type=float,
                    default=None,
                    required=True,
                    help='Number of hours for slurm job')

arser.add_argument('--out-dir',
                    type=str,
                    default=None,
                    required=True,
                    help='Output directory')

parser.add_argument('--drq',
                    type=str,
                    default=None,
                    required=True,
                    help='Catalog of objects in DRQ format')

parser.add_argument('--in-dir',
                    type=str,
                    default=None,
                    required=True,
                    help='Directory containing spectra files')

parser.add_argument('--zqso-min',
                    type=float,
                    default=None,
                    required=False,
                    help='Lower limit on quasar redshift from drq')

parser.add_argument('--lambda-min',
                    type=float,
                    default=3600.,
                    required=False,
                    help='Lower limit on observed wavelength [Angstrom]')

parser.add_argument('--lambda-max',
                    type=float,
                    default=5500.,
                    required=False,
                    help='Upper limit on observed wavelength [Angstrom]')

parser.add_argument('--lambda-rest-min',
                    type=float,
                    default=1040.,
                    required=False,
                    help='Lower limit on rest frame wavelength [Angstrom]')

parser.add_argument('--lambda-rest-max',
                    type=float,
                    default=1200.,
                    required=False,
                    help='Upper limit on rest frame wavelength [Angstrom]')

parser.add_argument('--rebin',
                    type=int,
                    default=3,
                    required=False,
                    help=('Rebin wavelength grid by combining this number '
                          'of adjacent pixels (ivar weight)'))

parser.add_argument('--npix-min',
                    type=int,
                    default=50,
                    required=False,
                    help='Minimum of rebined pixels')

parser.add_argument('--nproc',
                    type=int,
                    default=None,
                    required=False,
                    help='Number of processors')

args = parser.parse_args()

## Make the script text.
time = submit_utils.nh_to_hhmmss(args.slurm_hours)
text = submit_utils.make_header(queue='regular',nnodes=1,time=time,job_name='picca_deltas',err_file='run-picca-deltas-%j.err',out_file='run-picca-deltas-%j.out')
text += '\n\n'
text += 'export OMP_NUM_THREADS=1\n'
text += 'srun -n 1 -c 64 picca_deltas.py --in-dir {} --drq {} --out-dir {} --mode desi --iter-out-prefix {}/iter --log {}/picca_deltas.log --nproc {} --zqso-min {} --lambda-min {} --lambda-max {} --lambda-rest-min {} --lambda-rest-max {} --rebin {} --npix-min {}'.format(args.in_dir,args.drq,args.out_dir,args.out_dir,args.out_dir,args.nproc,args.zqso_min,args.lambda_min,args.lambda_max,args.lambda_rest_min,args.lambda_rest_max,args.rebin,args.npix_min)

## Write the script.
with open(args.slurm_script,'w') as f:
    f.write(text)
submit_utils.make_file_executable(args.slurm_script)

## Submit the job.
call('sbatch {}'.format(args.slurm_script))
