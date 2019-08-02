import numpy as np
from subprocess import call
import argparse
import os

from lyacolore import submit_utils

################################################################################

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--base-dir', type = str, default = None, required=False,
                    help = 'directory containing the master file')

parser.add_argument('--v-maj', type = int, default = 9, required=False,
                    help = 'major version of lyacolore realisations')

parser.add_argument('--v-min', type = int, default = 0, required=False,
                    help = 'minor version of lyacolore realisations')

parser.add_argument('--v-realisations', type = int, default = 0, required=False,
                    help = 'realisation numbers of lyacolore realisations', nargs='*')

parser.add_argument('--nproc', type = int, default = 32, required=False,
                    help = 'number of processes to use')

parser.add_argument('--nside', type = int, default = 16, required=False,
                    help = 'HEALPix nside for output files (must be 2^n)')

parser.add_argument('--fid-Om', type=float, default=0.315, required=False,
                    help='Omega_matter(z=0) of fiducial LambdaCDM cosmology')

parser.add_argument('--fid-Or', type=float, default=0., required=False,
                    help='Omega_radiation(z=0) of fiducial LambdaCDM cosmology')

args = parser.parse_args()

################################################################################

data_dir = args.base_dir + '/data/'
add_data_dir = data_dir + '/additional_data/'
submit_utils.check_dir(add_data_dir)
zbins = [(0.,10.)]

global_job_info = {'zbins':     zbins,
                   'options':   {'fid-Om':    args.fid_Om,
                                 'fid-Or':    args.fid_Or,
                                 'nside':     args.nside,
                                 'nproc':     args.nproc,
                                 },
                   }

for v_rea in args.v_realisations:

    for zbin in zbins:

        zmin = zbin[0]
        zmax = zbin[1]

        #Check that the directories are set up correctly.
        ver = 'v{}.{}.{}'.format(args.v_maj,args.v_min,v_rea)
        print('\nRunning analysis for version {}:'.format(ver))
        deltas_dir = data_dir + '/picca_input/' + ver + '/deltas_Lyb_metals/'
        submit_utils.check_dir(add_data_dir+'/run_files/')
        submit_utils.check_dir(add_data_dir+'/scripts/')

        #Determine the job parameters and names to make the header info.
        time = 20.
        job_name = 'make_metal_dmat_{}_{}_{}'.format(ver,zmin,zmax)
        err_file = 'make-metal-dmat-{}-{}-{}-%j.err'.format(ver,zmin,zmax)
        out_file = 'make-metal-dmat-{}-{}-{}-%j.out'.format(ver,zmin,zmax)

        #Construct the job info directory.
        header_info = {'queue':    'regular',
                       'time':     '{}:00:00'.format(int(time)),
                       'job_name': job_name,
                       'err_file': err_file,
                       'out_file': out_file,
                       }

        options = {'in-dir':    deltas_dir,
                   'out':       add_data_dir+'/metal_dmat_{}_{}_{}.fits.gz'.format(ver,zmin,zmax),
                   'rej':       0.999,
                   'abs-igm':   'SiII\(1260\) SiIII\(1207\) SiII\(1193\) SiII\(1190\)',
                   'z-cut-min': zmin,
                   'z-cut-max': zmax,
                   }

        metal_dmat_job_info = {'dir':            add_data_dir,
                               'header_info':    header_info,
                               'picca_script':   'picca_metal_dmat.py',
                               'options':        options,
                               'run_script':     'run_metal_dmat_{}_{}.sh'.format(zmin,zmax),
                               }

        #Submit the job.
        submit_utils.run_picca_job(metal_dmat_job_info,global_job_info['options'])
