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

parser.add_argument('--run-all', action="store_true", default = False, required=False,
                    help = 'run all correlations')

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
if args.run_all:
    args.run_lya_auto = True
    args.run_qso_auto = True
    args.run_dla_auto = True
    args.run_lya_aa_auto = True
    args.run_lya_qso_cross = True
    args.run_lya_dla_cross = True
    args.run_qso_dla_cross = True

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

#Function to check that a directory exists.
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

#Function to add an option to a command.
def add_to_command(command,extra,ns=False):
    if ns:
        return command + extra
    else:
        return command + ' ' + extra

#Function to interpret the job information and global options, and submit a job
#suitably to the queue.
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

    command = 'command="'
    command = add_to_command(command,job_info['picca_script'],ns=True)
    options = job_info['options']
    for key in options:
        command = add_to_command(command,'--{} {}'.format(key,options[key]))
    for key in global_options:
        command = add_to_command(command,'--{} {}'.format(key,global_options[key]))
    command = add_to_command(command,'"')
    command = add_to_command(command,'\nsrun -N 1 -n 1 -c {} $command\n'.format(args.nproc))

    #Make the run script.
    run_script_text = header + command
    run_script_path = dir+'/scripts/'+job_info['run_script']
    run_script = open(run_script_path,'w+')
    run_script.write(run_script_text)
    run_script.close()
    print(' -> -> -> job script written: {}'.format(job_info['run_script']))

    #Send the run script.
    print(' -> -> -> sending job to queue...')
    retcode = call('sbatch {}'.format(run_script_path),shell=True)
    print(' ')

    return

#Functions to construct the job info dictionaries.
def make_lya_auto_job_info(meas_dir,ver,zmin,zmax,deltas_dir):

    dir = meas_dir + '/lya_auto/'
    out_dir = dir + '/correlations/'
    check_dir(dir)

    header_info = {'queue':    'regular',
                   'time':     '12:00:00',
                   'job_name': 'run_lya_auto_{}_{}_{}'.format(ver,zmin,zmax),
                   'err_file': 'lya_auto_{}_{}_{}_%j.err'.format(ver,zmin,zmax),
                   'out_file': 'lya_auto_{}_{}_{}_%j.out'.format(ver,zmin,zmax),
                   }

    options = {'in-dir':        deltas_dir,
               'out':           out_dir+'/cf_lya_auto_{}_{}.fits.gz'.format(zmin,zmax),
               'no-project':    '',
               'z-cut-min':     zmin,
               'z-cut-max':     zmax,
               }

    lya_auto_job_info = {'dir':             dir,
                         'header_info':     header_info,
                         'picca_script':    'picca_cf.py',
                         'options':         options,
                         'run_script':      'run_lya_auto_{}_{}.sh'.format(zmin,zmax),
                         }

    return lya_auto_job_info

def make_qso_auto_job_info(meas_dir,ver,zmin,zmax,zcat,zcat_rand,corr_type='DD'):

    dir = meas_dir + '/qso_auto/'
    out_dir = dir + '/correlations/'
    check_dir(dir)

    if corr_type == 'DD':
        time = '01:30:00'
    elif corr_type == 'DR' or corr_type == 'RD':
        time = '03:00:00'
    elif corr_type == 'RR':
        time = '06:00:00'

    header_info = {'queue':    'regular',
                   'time':     time,
                   'job_name': 'run_qso_auto_{}_{}_{}_{}'.format(corr_type,ver,zmin,zmax),
                   'err_file': 'qso_auto_{}_{}_{}_{}_%j.err'.format(corr_type,ver,zmin,zmax),
                   'out_file': 'qso_auto_{}_{}_{}_{}_%j.out'.format(corr_type,ver,zmin,zmax),
                   }

    if corr_type == 'DD':
        options = {'drq':           zcat,
                   'z-evol-obj':    1.44,
                   }
    elif corr_type == 'DR':
        options = {'drq':           zcat,
                   'drq2':          zcat_rand,
                   'z-evol-obj':    1.44,
                   'z-evol-obj2':   1.44,
                   }
    elif corr_type == 'RD':
        options = {'drq':           zcat_rand,
                   'drq2':          zcat,
                   'z-evol-obj':    1.44,
                   'z-evol-obj2':   1.44,
                   }
    elif corr_type == 'RR':
        options = {'drq':           zcat_rand,
                   'z-evol-obj':    1.44,
                   }

    options = {**options,
               'out':           out_dir+'co_qso_auto_{}_{}_{}.fits.gz'.format(corr_type,zmin,zmax),
               'type-corr':     corr_type,
               'z-cut-min':     zmin,
               'z-cut-max':     zmax,
               }

    qso_auto_job_info = {'dir':             dir,
                         'header_info':     header_info,
                         'picca_script':    'picca_co.py',
                         'options':         options,
                         'run_script':      'run_qso_auto_{}_{}_{}.sh'.format(corr_type,zmin,zmax),
                         }

    return qso_auto_job_info

def make_dla_auto_job_info(meas_dir,ver,zmin,zmax,zcat,zcat_rand,corr_type='DD'):

    dir = meas_dir + '/dla_auto/'
    out_dir = dir + '/correlations/'
    check_dir(dir)

    if corr_type == 'DD':
        time = '01:00:00'
    elif corr_type == 'DR' or corr_type == 'RD':
        time = '01:30:00'
    elif corr_type == 'RR':
        time = '02:30:00'

    header_info = {'queue':    'regular',
                   'time':     time,
                   'job_name': 'run_dla_auto_{}_{}_{}_{}'.format(corr_type,ver,zmin,zmax),
                   'err_file': 'dla_auto_{}_{}_{}_{}_%j.err'.format(corr_type,ver,zmin,zmax),
                   'out_file': 'dla_auto_{}_{}_{}_{}_%j.out'.format(corr_type,ver,zmin,zmax),
                   }

    if corr_type == 'DD':
        options = {'drq':           zcat,
                   'z-evol-obj':    0.0,
                   }
    elif corr_type == 'DR':
        options = {'drq':           zcat,
                   'drq2':          zcat_rand,
                   'z-evol-obj':    0.0,
                   'z-evol-obj2':   0.0,
                   }
    elif corr_type == 'RD':
        options = {'drq':           zcat_rand,
                   'drq2':          zcat,
                   'z-evol-obj':    0.0,
                   'z-evol-obj2':   0.0,
                   }
    elif corr_type == 'RR':
        options = {'drq':           zcat_rand,
                   'z-evol-obj':    0.0,
                   }

    options = {**options,
               'out':           out_dir+'/co_dla_auto_{}_{}_{}.fits.gz'.format(corr_type,zmin,zmax),
               'type-corr':     corr_type,
               'z-cut-min':     zmin,
               'z-cut-max':     zmax,
               }

    dla_auto_job_info = {'dir':             dir,
                         'header_info':     header_info,
                         'picca_script':    'picca_co.py',
                         'options':         options,
                         'run_script':      'run_dla_auto_{}_{}_{}.sh'.format(corr_type,zmin,zmax),
                         }

    return dla_auto_job_info

def make_lya_aa_auto_job_info(meas_dir,ver,zmin,zmax,deltas_dir):

    dir = meas_dir + '/lya_aa_auto/'
    out_dir = dir + '/correlations/'
    check_dir(dir)

    header_info = {'queue':    'regular',
                   'time':     '12:00:00',
                   'job_name': 'run_lya_aa_auto_{}_{}_{}'.format(ver,zmin,zmax),
                   'err_file': 'lya_aa_auto_{}_{}_{}_%j.err'.format(ver,zmin,zmax),
                   'out_file': 'lya_aa_auto_{}_{}_{}_%j.out'.format(ver,zmin,zmax),
                   }

    options = {'in-dir':        deltas_dir,
               'out':           out_dir+'cf_lya_aa_auto_{}_{}.fits.gz'.format(zmin,zmax),
               'no-project':    '',
               'z-cut-min':     zmin,
               'z-cut-max':     zmax,
               }

    lya_aa_auto_job_info = {'dir':             dir,
                         'header_info':     header_info,
                         'picca_script':    'picca_cf.py',
                         'options':         options,
                         'run_script':      'run_lya_aa_auto_{}_{}.sh'.format(zmin,zmax),
                         }

    return lya_aa_auto_job_info

def make_lya_qso_cross_job_info(meas_dir,ver,zmin,zmax,deltas_dir,zcat,zcat_rand,cat_type='D'):

    dir = meas_dir + '/lya_qso_cross/'
    out_dir = dir + '/correlations/'
    check_dir(dir)

    if cat_type == 'D':
        time = '04:00:00'
    elif cat_type == 'R':
        time = '08:00:00'

    header_info = {'queue':    'regular',
                   'time':     time,
                   'job_name': 'run_lya_qso_cross_{}_{}_{}_{}'.format(cat_type,ver,zmin,zmax),
                   'err_file': 'lya_qso_cross_{}_{}_{}_{}_%j.err'.format(cat_type,ver,zmin,zmax),
                   'out_file': 'lya_qso_cross_{}_{}_{}_{}_%j.out'.format(cat_type,ver,zmin,zmax),
                   }

    if cat_type == 'D':
        options = {'drq':           zcat,
                   }
    elif cat_type == 'R':
        options = {'drq':           zcat_rand,
                   }

    options = {**options,
               'in-dir':                    deltas_dir,
               'out':                       out_dir+'xcf_lya_qso_cross_{}_{}_{}.fits.gz'.format(cat_type,zmin,zmax),
               'z-evol-obj':                1.44,
               'z-cut-min':                 zmin,
               'z-cut-max':                 zmax,
               'no-project':                '',
               'no-remove-mean-lambda-obs': '',
               }

    lya_qso_cross_job_info = {'dir':            dir,
                              'header_info':    header_info,
                              'picca_script':   'picca_xcf.py',
                              'options':        options,
                              'run_script':     'run_lya_qso_cross_{}_{}_{}.sh'.format(cat_type,zmin,zmax),
                              }

    return lya_qso_cross_job_info

def make_lya_dla_cross_job_info(meas_dir,ver,zmin,zmax,deltas_dir,zcat,zcat_rand,cat_type='D'):

    dir = meas_dir + '/lya_dla_cross/'
    out_dir = dir + '/correlations/'
    check_dir(dir)

    if cat_type == 'D':
        time = '04:00:00'
    elif cat_type == 'R':
        time = '08:00:00'

    header_info = {'queue':    'regular',
                   'time':     time,
                   'job_name': 'run_lya_dla_cross_{}_{}_{}_{}'.format(cat_type,ver,zmin,zmax),
                   'err_file': 'lya_dla_cross_{}_{}_{}_{}_%j.err'.format(cat_type,ver,zmin,zmax),
                   'out_file': 'lya_dla_cross_{}_{}_{}_{}_%j.out'.format(cat_type,ver,zmin,zmax),
                   }

    if cat_type == 'D':
        options = {'drq':   zcat,
                   }
    elif cat_type == 'R':
        options = {'drq':   zcat_rand,
                   }

    options = {**options,
               'in-dir':                    deltas_dir,
               'out':                       out_dir+'xcf_lya_dla_cross_{}_{}_{}.fits.gz'.format(cat_type,zmin,zmax),
               'z-evol-obj':                0.0,
               'z-cut-min':                 zmin,
               'z-cut-max':                 zmax,
               'no-project':                '',
               'no-remove-mean-lambda-obs': '',
               }

    lya_dla_cross_job_info = {'dir':            dir,
                              'header_info':    header_info,
                              'picca_script':   'picca_xcf.py',
                              'options':        options,
                              'run_script':     'run_lya_dla_cross_{}_{}_{}.sh'.format(cat_type,zmin,zmax),
                              }

    return lya_dla_cross_job_info

def make_qso_dla_cross_job_info(meas_dir,ver,zmin,zmax,zcat_qso,zcat_qso_rand,zcat_dla,zcat_dla_rand,corr_type='DD'):

    dir = meas_dir + '/qso_dla_cross/'
    out_dir = dir + '/correlations/'
    check_dir(dir)

    if corr_type == 'DD':
        time = '01:30:00'
    elif corr_type == 'DR' or corr_type == 'RD':
        time = '03:00:00'
    elif corr_type == 'RR':
        time = '06:00:00'

    header_info = {'queue':    'regular',
                   'time':     time,
                   'job_name': 'run_qso_dla_cross_{}_{}_{}_{}'.format(corr_type,ver,zmin,zmax),
                   'err_file': 'qso_dla_cross_{}_{}_{}_{}_%j.err'.format(corr_type,ver,zmin,zmax),
                   'out_file': 'qso_dla_cross_{}_{}_{}_{}_%j.out'.format(corr_type,ver,zmin,zmax),
                   }

    if corr_type == 'DD':
        options = {'drq':           zcat_qso,
                   'drq2':          zcat_dla,
                   'z-evol-obj':    1.44,
                   'z-evol-obj2':   0.0,
                   }
    elif corr_type == 'DR':
        options = {'drq':           zcat_qso,
                   'drq2':          zcat_dla_rand,
                   'z-evol-obj':    1.44,
                   'z-evol-obj2':   0.0,
                   }
    elif corr_type == 'RD':
        options = {'drq':           zcat_qso_rand,
                   'drq2':          zcat_dla,
                   'z-evol-obj':    1.44,
                   'z-evol-obj2':   0.0,
                   }
    elif corr_type == 'RR':
        options = {'drq':           zcat_qso_rand,
                   'drq2':          zcat_dla_rand,
                   'z-evol-obj':    1.44,
                   'z-evol-obj2':   0.0,
                   }

    options = {**options,
               'out':           out_dir+'co_qso_dla_cross_{}_{}_{}.fits.gz'.format(corr_type,zmin,zmax),
               'type-corr':     corr_type,
               'z-cut-min':     zmin,
               'z-cut-max':     zmax,
               }

    qso_dla_cross_job_info = {'dir':             dir,
                              'header_info':     header_info,
                              'picca_script':    'picca_co.py',
                              'options':         options,
                              'run_script':      'run_qso_dla_cross_{}_{}_{}.sh'.format(corr_type,zmin,zmax),
                              }

    return qso_dla_cross_job_info

################################################################################

global_job_info = {'zbins': [(0.0,2.2),(2.2,2.6),(2.6,3.0),(3.0,10.0)],
                   'options': {'fid-Om':    args.fid_Om,
                               'fid-Or':    args.fid_Or,
                               'nside':     args.nside,
                               'nproc':     args.nproc,
                               },
                   }

njobs = 0
for v_rea in args.v_realisations:

    ver = 'v{}.{}.{}'.format(args.v_maj,args.v_min,v_rea)
    print('\nRunning analysis for version {}:'.format(ver))

    #Check that the directories are constructed properly.
    ac_dir = a_dir+'/correlation_functions/'
    check_dir(ac_dir)
    acv_dir = ac_dir+'/'+ver+'/'
    check_dir(acv_dir)
    acvm_dir = acv_dir+'/measurements/'
    check_dir(acvm_dir)

    #Define the location variables for this version.
    lya_deltas_loc = '{}/data/picca_input/{}/deltas_0.5/'.format(args.base_dir,ver)
    lya_aa_deltas_loc = '{}/data/picca_input/{}/deltas_0.5_Lyb_metals/'.format(args.base_dir,ver)
    zcat_qso_loc = '{}/data/picca_input/{}/zcat_0.5.fits'.format(args.base_dir,ver)
    zcat_dla_loc = '{}/data/picca_input/{}/zcat_DLA_0.5.fits'.format(args.base_dir,ver)
    zcat_qso_rand_loc = '{}/data/picca_input/{}/zcat_0.1_randoms.fits'.format(args.base_dir,ver)
    zcat_dla_rand_loc = '{}/data/picca_input/{}/zcat_DLA_0.1_randoms.fits'.format(args.base_dir,ver)

    for zbin in global_job_info['zbins']:

        zmin = zbin[0]
        zmax = zbin[1]

        if args.run_lya_auto:

            print(' -> setting up lya auto correlation:')
            lya_auto_job_info = make_lya_auto_job_info(acvm_dir,ver,zmin,zmax,lya_deltas_loc)
            run_picca_job(lya_auto_job_info,global_job_info['options'])
            njobs += 1

        if args.run_qso_auto:

            print(' -> setting up qso auto correlation:')
            for corr_type in ['DD','RD','DR','RR']:
                print(' -> -> configuring {} run'.format(corr_type))
                qso_auto_job_info = make_qso_auto_job_info(acvm_dir,ver,zmin,zmax,zcat_qso_loc,zcat_qso_rand_loc,corr_type=corr_type)
                run_picca_job(qso_auto_job_info,global_job_info['options'])
                njobs += 1

        if args.run_dla_auto:

            print(' -> setting up dla auto correlation:')
            for corr_type in ['DD','RD','DR','RR']:
                print(' -> -> configuring {} run'.format(corr_type))
                dla_auto_job_info = make_dla_auto_job_info(acvm_dir,ver,zmin,zmax,zcat_dla_loc,zcat_dla_rand_loc,corr_type=corr_type)
                run_picca_job(dla_auto_job_info,global_job_info['options'])
                njobs += 1

        if args.run_lya_aa_auto:

            print(' -> setting up lya + all absorbers auto correlation:')
            lya_aa_auto_job_info = make_lya_aa_auto_job_info(acvm_dir,ver,zmin,zmax,lya_aa_deltas_loc)
            run_picca_job(lya_aa_auto_job_info,global_job_info['options'])
            njobs += 1

        if args.run_lya_qso_cross:

            print(' -> setting up lya qso cross correlation:')
            for cat_type in ['D','R']:
                print(' -> -> configuring run using {} catalog'.format(cat_type))
                lya_qso_cross_job_info = make_lya_qso_cross_job_info(acvm_dir,ver,zmin,zmax,lya_deltas_loc,zcat_qso_loc,zcat_qso_rand_loc,cat_type=cat_type)
                run_picca_job(lya_qso_cross_job_info,global_job_info['options'])
                njobs += 1

        if args.run_lya_dla_cross:

            print(' -> setting up lya dla cross correlation:')
            for cat_type in ['D','R']:
                print(' -> -> configuring run using {} catalog'.format(cat_type))
                lya_dla_cross_job_info = make_lya_dla_cross_job_info(acvm_dir,ver,zmin,zmax,lya_deltas_loc,zcat_dla_loc,zcat_dla_rand_loc,cat_type=cat_type)
                run_picca_job(lya_dla_cross_job_info,global_job_info['options'])
                njobs += 1

        if args.run_qso_dla_cross:

            print(' -> setting up qso dla cross correlation:')
            for corr_type in ['DD','RD','DR','RR']:
                print(' -> -> configuring {} run'.format(corr_type))
                qso_dla_cross_job_info = make_qso_dla_cross_job_info(acvm_dir,ver,zmin,zmax,zcat_qso_loc,zcat_qso_rand_loc,zcat_dla_loc,zcat_dla_rand_loc,corr_type=corr_type)
                run_picca_job(qso_dla_cross_job_info,global_job_info['options'])
                njobs += 1

################################################################################

print('\nAll analyses for all realisations sent to the queue (total {} jobs).'.format(njobs))

################################################################################
