import numpy as np
from subprocess import call
import argparse
import os
import fitsio

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
def run_picca_job(job_info,global_options,no_submit=False):

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

    if not no_submit:
        #Send the run script.
        retcode = call('sbatch {}'.format(run_script_path),shell=True)

    return

#Function to stack correlation files by simply concatenating the lists of
#subsamples in each correltion.
def concatenate_subsamples(fi,fout,corr_type):

    #Set up the data structures to store the information.
    cor_data = []
    cor = {}
    h = fitsio.FITS(fi[0])
    attri = {}
    for k in ['RP','RT','Z','NB']:
        attri[k] = np.zeros(h[1][k][:].shape)

    #Assume that the headers in the different files are all the same (and correct)
    head = h[1].read_header()
    head2 = h[2].read_header()
    h.close()

    #Loop through the files, adding the data to our overarching structures.
    for f in fi:
        h = fitsio.FITS(f)
        for k in ['RP','RT','Z']:
            attri[k] += h[1][k][:] * h[1]['NB'][:]
        attri['NB'] += h[1]['NB'][:]
        cor_data += [h[2][:]]
        h.close()

    #Ensure that data are correctly normalised.
    for k in ['RP','RT','Z']:
        attri[k] /= attri['NB']

    if corr_type in ['cf','xcf']:
        for k in ['WE','DA']:
            cor[k] = np.concatenate([cd[k] for cd in cor_data],axis=0)
    elif corr_type in ['co']:
        for k in ['WE','NB']:
            cor[k] = np.concatenate([cd[k] for cd in cor_data],axis=0)

    #Give shift the HEALPID of the pixels in order to avoid duplication.
    shift = int(10**(np.floor(np.log10(np.max(np.concatenate([cd['HEALPID'] for cd in cor_data])))) + 1))
    cor['HEALPID'] = np.concatenate([cd['HEALPID']+i*shift for i,cd in enumerate(cor_data)],axis=0)

    #Construct the output file in the standard picca format.
    out = fitsio.FITS(fout,'rw',clobber=True)
    names = ['RP','RT','Z','NB']
    out.write([attri[k] for k in names],names=names,
        comment=['R-parallel','R-transverse','Redshift','Number of pairs'],
        units=['h^-1 Mpc','h^-1 Mpc','',''],
        header=head,extname='ATTRI')
    head2 = [{'name':'HLPXSCHM','value':'RING','comment':'Healpix scheme'},
             {'name':'HIDSHIFT','value':shift,'comment':'Shift unit applied to HEALPix IDs'}]
    if corr_type in ['cf','xcf']:
        names2 = ['HEALPID','WE','DA']
        comment2 = ['Healpix index', 'Sum of weight', 'Correlation']
    elif corr_type in ['co']:
        names2 = ['HEALPID','WE','NB']
        comment2 = ['Healpix index', 'Sum of weight', 'Number of pairs']
    out.write([cor[k] for k in names2],names=names2,
        comment=comment2,
        header=head2,extname='COR')

    out.close

    return

def nh_to_hhmmss(nh):
    hh = int(np.floor(nh))
    remh = nh - hh
    mm = int(np.floor(remh*60))
    remm = remh*60 - mm
    ss = int(np.ceil(remm*60))
    if ss==60:
        ss = 0
        mm += 1
    hhmmss = '{:02d}:{:02d}:{:02d}'.format(hh,mm,ss)
    return hhmmss

def make_permission_group_desi(dir):
    call(['chgrp', '-R', 'desi', dir])
    return

def make_file_executable(f):
    call(['chmod', 'ug+x', f])
    return

def make_analysis_dir(dirname):
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    make_permission_group_desi(dirname)

    for ext in ['data','correlations','fits']:
        extdirname = os.path.join(dirname,ext)
        if not os.path.isdir(extdirname):
            os.mkdir(extdirname)
        make_permission_group_desi(extdirname)

    deltasdirname = os.path.join(dirname,'data/deltas')
    if not os.path.isdir(deltasdirname):
        os.mkdir(deltasdirname)
    make_permission_group_desi(deltasdirname)

    return

class AnalysisDir:

    def __init__(self, dirloc, dirname, datadirname='data', corrdirname='correlations', fitsdirname='fits', deltadirname='deltas'):

        self.dirloc = dirloc
        self.dirname = dirname
        self.datadirname = datadirname
        self.corrdirname = corrdirname
        self.fitsdirname = fitsdirname
        self.deltadirname = deltadirname

        if not os.path.isdir(self.dirloc):
            os.mkdir(self.dirloc)
        make_permission_group_desi(self.dirloc)

        self.dir = os.path.join(self.dirloc,self.dirname)
        if not os.path.isdir(self.dir):
            os.mkdir(self.dir)
        make_permission_group_desi(self.dir)

        self.datadir = os.path.join(self.dir,self.datadirname)
        if not os.path.isdir(self.datadir):
            os.mkdir(self.datadir)
        make_permission_group_desi(self.datadir)

        self.corrdir = os.path.join(self.dir,self.corrdirname)
        if not os.path.isdir(self.corrdir):
            os.mkdir(self.corrdir)
        make_permission_group_desi(self.corrdir)

        self.fitsdir = os.path.join(self.dir,self.fitsdirname)
        if not os.path.isdir(self.fitsdir):
            os.mkdir(self.fitsdir)
        make_permission_group_desi(self.fitsdir)

        self.deltadir = os.path.join(self.datadir,self.deltadirname)
        if not os.path.isdir(self.deltadir):
            os.mkdir(self.deltadir)
        make_permission_group_desi(self.deltadir)

        return
