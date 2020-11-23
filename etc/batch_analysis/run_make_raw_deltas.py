#!/usr/bin/env python

import argparse
import os

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

parser.add_argument('--out-dir',
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
text += 'converters.desi_convert_transmission_to_delta_files({}+"/drq.fits", {}+"/deltas/",in_dir={}, lambda_min={}, lambda_max={}, lambda_min_rest_frame={}, lambda_max_rest_frame={}, delta_log_lambda=3.e-4, nproc={})'.format(args.out_dir,args.out_dir,args.in_dir,args.lambda_min,args.lambda_max,,args.lambda_min_rest_frame,,args.lambda_max_rest_frame,args.nproc)

## Write the python script.
with open(args.python_script,'rw') as f:
    f.write(text)

## Make the slurm script text.
time = submit_utils.nh_to_hhmmss(args.slurm_hours)
text = submit_utils.make_header(queue='regular',nnodes=1,time=time,job_name='picca_deltas',err_file='run-picca-deltas-%j.err',out_file='run-picca-deltas-%j.out')
text += '\n\n'
text += 'export OMP_NUM_THREADS=1\n'
text += 'srun -n 1 -c 64 {}'.format(args.python_script)

## Write the slurm script.
with open(args.slurm_script,'rw') as f:
    f.write(text)

## Submit the job.
call('sbatch {}'.format(args.slurm_script))
