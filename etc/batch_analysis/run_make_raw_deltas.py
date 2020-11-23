#!/usr/bin/env python

import argparse
import os

from subprocess import call
from lyacolore import submit_utils

parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )

parser.add_argument('--python-script',
                    type=str,
                    default=None,
                    required=True,
                    help='Output python script')

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

parser.add_argument('--slurm-queue',
                    type=str,
                    default=None,
                    required=True,
                    help='Slurm queue to use')

parser.add_argument('-o','--out-dir',
                    type=str,
                    default=None,
                    required=True,
                    help='Output directory (inside which deltas dir will be made)')

parser.add_argument('--drq',
                    type=str,
                    default=None,
                    required=True,
                    help='Catalog of objects in DRQ format')

parser.add_argument('-i','--in-dir',
                    type=str,
                    default=None,
                    required=True,
                    help='Directory of raw mocks')

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

parser.add_argument('--nproc',
                    type=int,
                    default=None,
                    required=False,
                    help='Number of processors')

args = parser.parse_args()

## Make the python script.
text = '#!/usr/bin/env python\n\n'
text += 'from picca import converters\n\n'
text += 'converters.desi_convert_transmission_to_delta_files('
text += '"{}", '.format(args.drq)
text += '"{}", '.format(os.path.join(args.out_dir,'deltas/'))
text += 'in_dir="{}", '.format(args.in_dir)
text += 'lambda_min={}, '.format(args.lambda_min)
text += 'lambda_max={}, '.format(args.lambda_max)
text += 'lambda_min_rest_frame={}, '.format(args.lambda_rest_min)
text += 'lambda_max_rest_frame={}, '.format(args.lambda_rest_max)
text += 'delta_log_lambda=3.e-4, '
text += 'nproc={})\n'.format(args.nproc)

## Write the python script.
with open(args.python_script,'w') as f:
    f.write(text)
submit_utils.make_file_executable(args.python_script)

## Make the slurm script text.
time = submit_utils.nh_to_hhmmss(args.slurm_hours)
err_file = os.path.join(args.out_dir,'run-picca-deltas-%j.err')
out_file = os.path.join(args.out_dir,'run-picca-deltas-%j.out')
text = submit_utils.make_header(queue=args.slurm_queue,nnodes=1,time=time,job_name='raw_picca_deltas',err_file=err_file,out_file=out_file)
text += 'export OMP_NUM_THREADS=1\n'
text += 'srun -n 1 -c 64 {}\n'.format(args.python_script)

## Write the slurm script.
with open(args.slurm_script,'w') as f:
    f.write(text)
submit_utils.make_file_executable(args.slurm_script)

## Submit the job.
call(['sbatch',args.slurm_script])
