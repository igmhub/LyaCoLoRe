import numpy as np
from subprocess import call
import argparse
import os

#Function to make the header at the top of each run script.
def make_header(queue='regular',nnodes=1,time='01:00:00',job_name='run_script',err_file='run-%j.err',out_file='run-%j.out'):

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

    dir = job_info['dir']

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
    command = add_to_command(command,'\nsrun -N 1 -n 1 -c {} $command\n'.format(global_options['nproc']))

    #Make the run script.
    run_script_text = header + command
    run_script_path = dir+'/scripts/'+job_info['run_script']
    run_script = open(run_script_path,'w+')
    run_script.write(run_script_text)
    run_script.close()
    print(' -> -> job script written: {}'.format(job_info['run_script']))

    #Send the run script.
    retcode = call('sbatch {}'.format(run_script_path),shell=True)

    return
