import numpy as np
from subprocess import call
import argparse
import os

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

parser.add_argument('--run-lya-auto', action="store_true", default = False, required=False,
                    help = 'run the lya auto correlation')

parser.add_argument('--run-qso-auto', action="store_true", default = False, required=False,
                    help = 'run the qso auto correlation')

parser.add_argument('--run-dla-auto', action="store_true", default = False, required=False,
                    help = 'run the dla auto correlation')

parser.add_argument('--run-lya-aa-auto', action="store_true", default = False, required=False,
                    help = 'run the lya all absorber auto correlation')

parser.add_argument('--run-lya-qso-cross', action="store_true", default = False, required=False,
                    help = 'run the lya-qso cross correlation')

parser.add_argument('--run-lya-dla-cross', action="store_true", default = False, required=False,
                    help = 'run the lya-dla cross correlation')

parser.add_argument('--run-qso-dla-cross', action="store_true", default = False, required=False,
                    help = 'run the qso-dla cross correlation')

parser.add_argument('--fid-Om', type=float, default=0.315, required=False,
                    help='Omega_matter(z=0) of fiducial LambdaCDM cosmology')

parser.add_argument('--fid-Or', type=float, default=0., required=False,
                    help='Omega_radiation(z=0) of fiducial LambdaCDM cosmology')

args = parser.parse_args()

################################################################################

a_dir = args.base_dir+'/analysis/'

################################################################################
#Function to make the header at the top of each run script.
def make_header(queue='regular',nnodes=1,time='00:01:00',job_name='run_script',err_file='run-%j.err',out_file='run-%j.err'):

    header = ''

    header += '#!/bin/bash -l\n\n'

    header += '#SBATCH --partition {}\n'.format(queue)
    header += '#SBATCH --nodes {}\n'.format(nnodes)
    header += '#SBATCH --time {}\n'.format(time)
    header += '#SBATCH --job-name {}\n'.format(job_name)
    header += '#SBATCH --error {}\n'.format(err_file)
    header += '#SBATCH --output {}\n'.format(out_file)
    header += '#SBATCH -C haswell\n'
    header += '#SBATCH -A desi\n\n'

    header += 'umask 0002\n'
    header += 'export OMP_NUM_THREADS=64\n\n'

    return header

#For each realisation, for each correlation desired, create a job script and
#send it to the queue.
njobs = 0
zbins = [(0.0,2.2),(2.2,2.6),(2.6,3.0),(3.0,10.0)]
for v_rea in args.v_realisations:

    ver = 'v{}.{}.{}'.format(args.v_maj,args.v_min,v_rea)
    print('\nRunning analysis for version {}:'.format(ver))
    avc_dir = a_dir+'/correlation_functions/'+ver+'/measurements/'

    #Define the location variables.
    lya_deltas_loc = '{}/data/picca_input/{}/deltas_0.5/ '.format(args.base_dir,ver)
    lya_aa_deltas_loc = '{}/data/picca_input/{}/deltas_0.5_Lyb_metals/ '.format(args.base_dir,ver)
    zcat_qso_loc = '{}/data/picca_input/{}/zcat_0.5.fits/ '.format(args.base_dir,ver)
    zcat_dla_loc = '{}/data/picca_input/{}/zcat_DLA_0.5.fits/ '.format(args.base_dir,ver)
    zcat_qso_rand_loc = '{}/data/picca_input/{}/zcat_0.1_randoms.fits/ '.format(args.base_dir,ver)
    zcat_dla_rand_loc = '{}/data/picca_input/{}/zcat_DLA_0.1_randoms.fits/ '.format(args.base_dir,ver)

    if args.run_lya_auto:

        print(' -> setting up the lya auto correlation...')

        for zbin in zbins:
            zmin = zbin[0]
            zmax = zbin[1]

            lya_auto_dir = avc_dir+'/lya_auto/'
            lya_auto_file = 'cf_lya_auto_{}_{}.fits.gz'.format(zmin,zmax)
            try:
                os.mkdir(lya_auto_dir)
            except FileExistsError:
                pass
            try:
                os.mkdir(lya_auto_dir+'/scripts/')
            except FileExistsError:
                pass
            try:
                os.mkdir(lya_auto_dir+'/correlations/')
            except FileExistsError:
                pass
            try:
                os.mkdir(lya_auto_dir+'/run_files/')
            except FileExistsError:
                pass

            #Make the header.
            queue = 'regular'
            time = '00:12:00'
            job_name = 'run_lya_auto_{}_{}_{}'.format(ver,zmin,zmax)
            err_file = lya_auto_dir+'/run_files/lya_auto_{}_{}_{}_%j.err'.format(ver,zmin,zmax)
            out_file = lya_auto_dir+'/run_files/lya_auto_{}_{}_{}_%j.out'.format(ver,zmin,zmax)
            header = make_header(queue=queue,time=time,job_name=job_name,err_file=err_file,out_file=out_file)

            #Make the command.
            command = ''
            command += 'command="picca_cf.py '
            command += '--in-dir {} '.format(lya_deltas_loc)
            command += '--out {}/correlations/{} '.format(lya_auto_dir,lya_auto_file)
            command += '--fid-Om {} '.format(args.fid_Om)
            command += '--fid-Or {} '.format(args.fid_Or)
            command += '--no-project '
            command += '--nside {} '.format(args.nside)
            command += '--nproc {} '.format(args.nproc)
            command += '--z-cut-min {} '.format(zmin)
            command += '--z-cut-max {} '.format(zmax)
            command += '"'
            command += '\n'
            command += 'srun -N 1 -n 1 -c {} $command'.format(args.nproc)
            command += '\n'

            #Make the run script.
            run_script_text = header + command
            run_script_path = '{}/scripts/run_lya_auto_{}_{}.sh'.format(lya_auto_dir,zmin,zmax)
            run_script = open(run_script_path,'w+')
            run_script.write(run_script_text)
            run_script.close()
            print(' -> -> job script written to:{}'.format(run_script_path))

            #Send the run script.
            print(' -> -> sending job to queue...')
            retcode = call('sbatch {}'.format(run_script_path),shell=True)
            njobs += 1
            print(' ')

    if args.run_qso_auto:

        print(' -> setting up the qso auto correlation...')
        zbins = [(0.0,2.2),(2.2,2.6),(2.6,3.0),(3.0,10.0)]

        for zbin in zbins:
            zmin = zbin[0]
            zmax = zbin[1]

            qso_auto_dir = avc_dir+'/qso_auto/'
            try:
                os.mkdir(qso_auto_dir)
            except FileExistsError:
                pass
            try:
                os.mkdir(qso_auto_dir+'/scripts/')
            except FileExistsError:
                pass
            try:
                os.mkdir(qso_auto_dir+'/correlations/')
            except FileExistsError:
                pass
            try:
                os.mkdir(qso_auto_dir+'/run_files/')
            except FileExistsError:
                pass

            for corr_type in ['DD','DR','RD','RR']:

                qso_auto_file = 'co_qso_auto_{}_{}_{}.fits.gz'.format(corr_type,zmin,zmax)

                #Make the header.
                queue = 'debug'
                time = '00:00:30'
                job_name = 'run_qso_auto_{}_{}_{}_{}'.format(ver,corr_type,zmin,zmax)
                err_file = qso_auto_dir+'/run_files/qso_auto_{}_{}_{}_{}_%j.err'.format(ver,corr_type,zmin,zmax)
                out_file = qso_auto_dir+'/run_files/qso_auto_{}_{}_{}_{}_%j.out'.format(ver,corr_type,zmin,zmax)
                header = make_header(queue=queue,time=time,job_name=job_name,err_file=err_file,out_file=out_file)

                #Make the command.
                command = ''
                command += 'command="picca_co.py '
                if corr_type in ['DD','DR']:
                    command += '--drq {} '.format(zcat_qso_loc)
                if corr_type in ['RR']:
                    command += '--drq {} '.format(zcat_qso_rand_loc)
                if corr_type in ['RD']:
                    command += '--drq2 {} '.format(zcat_qso_loc)
                if corr_type in ['DR','RR']:
                    command += '--drq2 {} '.format(zcat_qso_rand_loc)
                command += '--out {}/correlations/{} '.format(qso_auto_dir,qso_auto_file)
                command += '--fid-Om {} '.format(args.fid_Om)
                command += '--fid-Or {} '.format(args.fid_Or)
                command += '--no-project '
                command += '--nside {} '.format(args.nside)
                command += '--nproc {} '.format(args.nproc)
                command += '--z-cut-min {} '.format(zmin)
                command += '--z-cut-max {} '.format(zmax)
                command += '--z-evol-obj 1.44 '
                if corr_type in ['DR','RR']:
                    command += '--z-evol-obj2 1.44 '
                command += '--type-corr {} '.format(corr_type)
                command += '"'
                command += '\n'
                command += 'srun -N 1 -n 1 -c {} $command'.format(args.nproc)
                command += '\n'

                #Make the run script.
                run_script_text = header + command
                run_script_path = '{}/scripts/run_qso_auto_{}_{}_{}.sh'.format(qso_auto_dir,corr_type,zmin,zmax)
                run_script = open(run_script_path,'w+')
                run_script.write(run_script_text)
                run_script.close()
                print(' -> -> job script written to:{}'.format(run_script_path))

                #Send the run script.
                print(' -> -> sending job to queue...')
                retcode = call('sbatch {}'.format(run_script_path),shell=True)
                njobs += 1
                print(' ')

    if args.run_lya_aa_auto:

        print(' -> setting up the lya_aa auto correlation...')
        zbins = [(0.0,2.2),(2.2,2.6),(2.6,3.0),(3.0,10.0)]

        for zbin in zbins:
            zmin = zbin[0]
            zmax = zbin[1]

            lya_aa_auto_dir = avc_dir+'/lya_aa_auto/'
            lya_aa_auto_file = 'cf_lya_aa_auto_{}_{}.fits.gz'.format(zmin,zmax)
            try:
                os.mkdir(lya_aa_auto_dir)
            except FileExistsError:
                pass
            try:
                os.mkdir(lya_aa_auto_dir+'/scripts/')
            except FileExistsError:
                pass
            try:
                os.mkdir(lya_aa_auto_dir+'/correlations/')
            except FileExistsError:
                pass
            try:
                os.mkdir(lya_aa_auto_dir+'/run_files/')
            except FileExistsError:
                pass

            #Make the header.
            queue = 'regular'
            time = '00:12:00'
            job_name = 'run_lya_aa_auto_{}_{}_{}'.format(ver,zmin,zmax)
            err_file = lya_aa_auto_dir+'/run_files/lya_aa_auto_{}_{}_{}_%j.err'.format(ver,zmin,zmax)
            out_file = lya_aa_auto_dir+'/run_files/lya_aa_auto_{}_{}_{}_%j.out'.format(ver,zmin,zmax)
            header = make_header(queue=queue,time='00:12:00',job_name=job_name,err_file=err_file,out_file=out_file)

            #Make the command.
            command = ''
            command += 'command="picca_cf.py '
            command += '--in-dir {}/data/picca_input/{}/deltas_0.5_Lyb_metals/ '.format(args.base_dir,ver)
            command += '--out {}/correlations/{} '.format(lya_aa_auto_dir,lya_aa_auto_file)
            command += '--fid-Om {} '.format(args.fid_Om)
            command += '--fid-Or {} '.format(args.fid_Or)
            command += '--no-project '
            command += '--nside {} '.format(args.nside)
            command += '--nproc {} '.format(args.nproc)
            command += '--z-cut-min {} '.format(zmin)
            command += '--z-cut-max {} '.format(zmax)
            command += '"'
            command += '\n'
            command += 'srun -N 1 -n 1 -c {} $command'.format(args.nproc)
            command += '\n'

            #Make the run script.
            run_script_text = header + command
            run_script_path = '{}/scripts/run_lya_aa_auto_{}_{}.sh'.format(lya_aa_auto_dir,zmin,zmax)
            run_script = open(run_script_path,'w+')
            run_script.write(run_script_text)
            run_script.close()
            print(' -> -> job script written to:{}'.format(run_script_path))

            #Send the run script.
            print(' -> -> sending job to queue...')
            retcode = call('sbatch {}'.format(run_script_path),shell=True)
            njobs += 1
            print(' ')

print('\nAll analyses for all realisations sent to the queue (total {} jobs).'.format(njobs))

################################################################################
