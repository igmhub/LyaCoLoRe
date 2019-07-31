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

a_dir = args.base_dir + '/analysis/'
picca_data_dir = args.base_dir + '/data/picca_input'

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

def check_dir(dir):

    try:
        os.mkdir(dir)
    except FileExistsError:
        pass

    return

#Function to make sure that the directories inside the correlation directory
#are set up properly.
def check_corr_dir(corr_dir):

    check_dir(corr_dir)
    check_dir(corr_dir+'/scripts/')
    check_dir(corr_dir+'/correlations/')
    check_dir(corr_dir+'/run_files/')

    return

def add_to_command(command,extra):
    return command + ' ' + extra

def run_picca_job(job_info,global_options):

    #Check the directory is set up.
    dir = job_info['dir']
    check_corr_dir(dir)

    #Make the header.
    header_info = job_info['header_info']
    header = make_header(queue=header_info['queue'],
                         time=header_info['time'],
                         job_name=header_info['job_name'],
                         err_file=dir+'/run_files/'+header_info['err_file'],
                         out_file=dir+'/run_files/'+header_info['out_file']
                         )

    command = 'command='
    command = add_to_command(command,job_info.picca_script)
    for key in job_info.options:
        command = add_to_command(command,'--{} {}'.format(key,job_info[key]))
    for key in global_options:
        command = add_to_command(command,'--{} {}'.format(key,job_info[key]))

    #Make the run script.
    run_script_text = header + command
    run_script_path = dir+'/scripts/'+job_info['run_script']
    run_script = open(run_script_path,'w+')
    run_script.write(run_script_text)
    run_script.close()
    print(' -> -> job script written to:{}'.format(run_script_path))

    #Send the run script.
    print(' -> -> sending job to queue...')
    retcode = call('sbatch {}'.format(run_script_path),shell=True)
    njobs += 1
    print(' ')

    return

def make_lya_auto_job_info(meas_dir,ver,zmin,zmax):

    dir = meas_dir + '/lya_auto/'
    check_dir(dir)
    
    header_info = {'queue':    'debug',
                   'time':     '00:05:00',
                   'job_name': 'run_lya_auto_{}_{}_{}'.format(ver,zmin,zmax),
                   'err_file': 'lya_auto_{}_{}_{}_%j.err'.format(ver,zmin,zmax),
                   'out_file': 'lya_auto_{}_{}_{}_%j.out'.format(ver,zmin,zmax),
                   }

    options = {'in-dir':        '{}/{}/deltas_0.5/ '.format(picca_data_dir,ver),
               'out':           'cf_lya_auto_{}_{}.fits.gz'.format(zmin,zmax),
               'no-project':    '',
               }

    lya_auto_job_info = {'dir':             dir,
                         'header_info':     header_info,
                         'picca_script':    'picca_cf.py',
                         'options':         options,
                         'run_script':      'run_lya_auto_{}_{}.sh'.format(zmin,zmax),
                         }

    return lya_auto_job_info

################################################################################

global_job_info = {'zbins': [(0.0,2.2),(2.2,2.6),(2.6,3.0),(3.0,10.0)],
                   'options': {'fid-Om':    args.fid_Om,
                               'fid-Or':    args.fid_Or,
                               'nside':     args.nside,
                               'nproc':     args.nproc,
                               },
                   }

for v_rea in args.v_realisations:

    ver = 'v{}.{}.{}'.format(args.v_maj,args.v_min,v_rea)
    print('\nRunning analysis for version {}:'.format(ver))
    ac_dir = a_dir+'/correlation_functions/'
    check_dir(ac_dir)
    acv_dir = ac_dir+'/'+ver+'/'
    check_dir(acv_dir)
    acvm_dir = acv_dir+'/measurements/'
    check_dir(acvm_dir)

    for zbin in global_job_info['zbins']:

        zmin = zbin[0]
        zmax = zbin[1]

        if args.run_lya_auto:

            lya_auto_job_info = make_lya_auto_job_info(acvm_dir,ver,zmin,zmax)
            run_picca_job(lya_auto_job_info,global_job_info['options'])




"""
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
    zcat_qso_loc = '{}/data/picca_input/{}/zcat_0.5.fits '.format(args.base_dir,ver)
    zcat_dla_loc = '{}/data/picca_input/{}/zcat_DLA_0.5.fits '.format(args.base_dir,ver)
    zcat_qso_rand_loc = '{}/data/picca_input/{}/zcat_0.1_randoms.fits '.format(args.base_dir,ver)
    zcat_dla_rand_loc = '{}/data/picca_input/{}/zcat_DLA_0.1_randoms.fits '.format(args.base_dir,ver)

    if args.run_lya_auto:

        print(' -> setting up the lya auto correlation...')

        for zbin in zbins:
            zmin = zbin[0]
            zmax = zbin[1]

            lya_auto_dir = avc_dir+'/lya_auto/'
            check_corr_dir(lya_auto_dir)
            lya_auto_file = 'cf_lya_auto_{}_{}.fits.gz'.format(zmin,zmax)

            #Make the header.
            queue = 'regular'
            time = '12:00:00'
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

        for zbin in zbins:
            zmin = zbin[0]
            zmax = zbin[1]

            qso_auto_dir = avc_dir+'/qso_auto/'
            check_corr_dir(qso_auto_dir)

            for corr_type in ['DD','DR','RD','RR']:

                qso_auto_file = 'co_qso_auto_{}_{}_{}.fits.gz'.format(corr_type,zmin,zmax)

                #Make the header.
                queue = 'regular'
                time = '04:00:00'
                job_name = 'run_qso_auto_{}_{}_{}_{}'.format(ver,corr_type,zmin,zmax)
                err_file = qso_auto_dir+'/run_files/qso_auto_{}_{}_{}_{}_%j.err'.format(ver,corr_type,zmin,zmax)
                out_file = qso_auto_dir+'/run_files/qso_auto_{}_{}_{}_{}_%j.out'.format(ver,corr_type,zmin,zmax)
                header = make_header(queue=queue,time=time,job_name=job_name,err_file=err_file,out_file=out_file)

                #Make the command.
                command = ''
                command += 'command="picca_co.py '
                if corr_type in ['DD','DR']:
                    command += '--drq {} '.format(zcat_qso_loc)
                if corr_type in ['RD','RR']:
                    command += '--drq {} '.format(zcat_qso_rand_loc)
                if corr_type in ['RD']:
                    command += '--drq2 {} '.format(zcat_qso_loc)
                if corr_type in ['DR']:
                    command += '--drq2 {} '.format(zcat_qso_rand_loc)
                command += '--out {}/correlations/{} '.format(qso_auto_dir,qso_auto_file)
                command += '--fid-Om {} '.format(args.fid_Om)
                command += '--fid-Or {} '.format(args.fid_Or)
                command += '--nside {} '.format(args.nside)
                command += '--nproc {} '.format(args.nproc)
                command += '--z-cut-min {} '.format(zmin)
                command += '--z-cut-max {} '.format(zmax)
                command += '--z-evol-obj 1.44 '
                if corr_type in ['DR','RD']:
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

    if args.run_dla_auto:

        print(' -> setting up the dla auto correlation...')

        for zbin in zbins:
            zmin = zbin[0]
            zmax = zbin[1]

            dla_auto_dir = avc_dir+'/dla_auto/'
            check_corr_dir(dla_auto_dir)

            for corr_type in ['DD','DR','RD','RR']:

                dla_auto_file = 'co_dla_auto_{}_{}_{}.fits.gz'.format(corr_type,zmin,zmax)

                #Make the header.
                queue = 'regular'
                time = '04:00:00'
                job_name = 'run_dla_auto_{}_{}_{}_{}'.format(ver,corr_type,zmin,zmax)
                err_file = dla_auto_dir+'/run_files/dla_auto_{}_{}_{}_{}_%j.err'.format(ver,corr_type,zmin,zmax)
                out_file = dla_auto_dir+'/run_files/dla_auto_{}_{}_{}_{}_%j.out'.format(ver,corr_type,zmin,zmax)
                header = make_header(queue=queue,time=time,job_name=job_name,err_file=err_file,out_file=out_file)

                #Make the command.
                command = ''
                command += 'command="picca_co.py '
                if corr_type in ['DD','DR']:
                    command += '--drq {} '.format(zcat_dla_loc)
                if corr_type in ['RD','RR']:
                    command += '--drq {} '.format(zcat_dla_rand_loc)
                if corr_type in ['RD']:
                    command += '--drq2 {} '.format(zcat_dla_loc)
                if corr_type in ['DR']:
                    command += '--drq2 {} '.format(zcat_dla_rand_loc)
                command += '--out {}/correlations/{} '.format(dla_auto_dir,dla_auto_file)
                command += '--fid-Om {} '.format(args.fid_Om)
                command += '--fid-Or {} '.format(args.fid_Or)
                command += '--nside {} '.format(args.nside)
                command += '--nproc {} '.format(args.nproc)
                command += '--z-cut-min {} '.format(zmin)
                command += '--z-cut-max {} '.format(zmax)
                command += '--z-evol-obj 0.0 '
                if corr_type in ['DR','RD']:
                    command += '--z-evol-obj2 0.0 '
                command += '--type-corr {} '.format(corr_type)
                command += '"'
                command += '\n'
                command += 'srun -N 1 -n 1 -c {} $command'.format(args.nproc)
                command += '\n'

                #Make the run script.
                run_script_text = header + command
                run_script_path = '{}/scripts/run_dla_auto_{}_{}_{}.sh'.format(dla_auto_dir,corr_type,zmin,zmax)
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
        print('Not set up yet!')

    if args.run_lya_qso_cross:

        print(' -> setting up the lya qso cross correlation...')

        for zbin in zbins:
            zmin = zbin[0]
            zmax = zbin[1]

            lya_qso_cross_dir = avc_dir+'/lya_qso_cross/'
            check_corr_dir(lya_qso_cross_dir)

            for obj_type in ['D','R']:

                lya_qso_cross_file = 'xcf_lya_qso_cross_{}_{}_{}.fits.gz'.format(obj_type,zmin,zmax)

                #Make the header.
                queue = 'regular'
                time = '4:00:00'
                job_name = 'run_lya_qso_cross_{}_{}_{}'.format(ver,zmin,zmax)
                err_file = lya_qso_cross_dir+'/run_files/lya_qso_cross_{}_{}_{}_{}_%j.err'.format(ver,obj_type,zmin,zmax)
                out_file = lya_qso_cross_dir+'/run_files/lya_qso_cross_{}_{}_{}_{}_%j.out'.format(ver,obj_type,zmin,zmax)
                header = make_header(queue=queue,time=time,job_name=job_name,err_file=err_file,out_file=out_file)

                #Make the command.
                command = ''
                command += 'command="picca_xcf.py '
                command += '--in-dir {} '.format(lya_deltas_loc)
                if obj_type == 'D':
                    command += '--drq {} '.format(zcat_qso_loc)
                elif obj_type == 'R':
                    command += '--drq {} '.format(zcat_qso_rand_loc)
                command += '--out {}/correlations/{} '.format(lya_qso_cross_dir,lya_qso_cross_file)
                command += '--fid-Om {} '.format(args.fid_Om)
                command += '--fid-Or {} '.format(args.fid_Or)
                command += '--no-project '
                command += '--no-remove-mean-lambda-obs '
                command += '--nside {} '.format(args.nside)
                command += '--nproc {} '.format(args.nproc)
                command += '--z-cut-min {} '.format(zmin)
                command += '--z-cut-max {} '.format(zmax)
                command += '--z-evol-obj 1.44 '
                command += '"'
                command += '\n'
                command += 'srun -N 1 -n 1 -c {} $command'.format(args.nproc)
                command += '\n'

                #Make the run script.
                run_script_text = header + command
                run_script_path = '{}/scripts/run_lya_qso_cross_{}_{}.sh'.format(lya_qso_cross_dir,zmin,zmax)
                run_script = open(run_script_path,'w+')
                run_script.write(run_script_text)
                run_script.close()
                print(' -> -> job script written to:{}'.format(run_script_path))

                #Send the run script.
                print(' -> -> sending job to queue...')
                retcode = call('sbatch {}'.format(run_script_path),shell=True)
                njobs += 1
                print(' ')

    if args.run_lya_dla_cross:
        print('Not set up yet!')

    if args.run_qso_dla_cross:
        print('Not set up yet!')
"""

print('\nAll analyses for all realisations sent to the queue (total {} jobs).'.format(njobs))

################################################################################
