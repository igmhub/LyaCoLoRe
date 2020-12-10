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
                    default='run_picca_deltas.sl',
                    required=False,
                    help='Output slurm script')

parser.add_argument('--slurm-hours',
                    type=float,
                    default=1.0,
                    required=False,
                    help='Number of hours for slurm job')

parser.add_argument('--slurm-queue',
                    type=str,
                    default='regular',
                    required=False,
                    help='Slurm queue to use')

parser.add_argument('-o','--out-dir',
                    type=str,
                    default=None,
                    required=True,
                    help='Path to output directory')

parser.add_argument('--deltas-dir',
                    type=str,
                    default=None,
                    required=True,
                    help='Name of directory to put deltas in')

parser.add_argument('--region-name',
                    type=str,
                    default=None,
                    required=True,
                    help='Name of region computing deltas in (lya/lyb)')

parser.add_argument('--drq',
                    type=str,
                    default=None,
                    required=True,
                    help='Catalog of objects in DRQ format')

parser.add_argument('-i','--in-dir',
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

submit_utils.check_dir(os.path.join(args.out_dir,'run_files'))
submit_utils.check_dir(os.path.join(args.out_dir,'scripts'))
args.slurm_script = os.path.join(args.out_dir,'scripts',args.slurm_script)

## Make the script text.
time = submit_utils.nh_to_hhmmss(args.slurm_hours)
err_file = os.path.join(args.out_dir,'run_files','run-picca-deltas-{}region-%j.err'.format(args.region_name))
out_file = os.path.join(args.out_dir,'run_files','run-picca-deltas-{}region-%j.out'.format(args.region_name))
text = submit_utils.make_header(queue=args.slurm_queue,nnodes=1,time=time,job_name='picca_deltas_{}region',err_file=err_file,out_file=out_file)
text += 'export OMP_NUM_THREADS=1\n'
text += 'srun -n 1 -c 64 picca_deltas.py '
text += '--in-dir {} '.format(args.in_dir)
text += '--drq {} '.format(args.drq)
text += '--out-dir {} '.format(os.path.join(args.out_dir,args.deltas_dir))
text += '--mode desi '
text += '--iter-out-prefix {}/iter_{}region '.format(args.out_dir,args.region_name)
text += '--log {}/picca_deltas_{}region.log '.format(args.out_dir,args.region_name)
if args.nproc is not None:
    text += '--nproc {} '.format(args.nproc)
if args.zqso_min is not None:
    text += '--zqso-min {} '.format(args.zqso_min)
text += '--lambda-min {} '.format(args.lambda_min)
text += '--lambda-max {} '.format(args.lambda_max)
text += '--lambda-rest-min {} '.format(args.lambda_rest_min)
text += '--lambda-rest-max {} '.format(args.lambda_rest_max)
text += '--rebin {} '.format(args.rebin)
text += '--npix-min {}\n'.format(args.npix_min)

## Write the script.
with open(args.slurm_script,'w') as f:
    f.write(text)
submit_utils.make_file_executable(args.slurm_script)

## Submit the job.
call(['sbatch',args.slurm_script])
