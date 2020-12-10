#!/usr/bin/env python

import numpy as np
from subprocess import call
import argparse
import os

from lyacolore import submit_utils

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--deltas-dir-lya', type = str, default = None, required=False,
                    help = 'directory containing the deltas for lya region')

parser.add_argument('--deltas-dir-lyb', type = str, default = None, required=False,
                    help = 'directory containing the deltas for lyb region')

parser.add_argument('--drq-qso', type = str, default = None, required=False,
                    help = 'qso catalogue')

parser.add_argument('--drq-dla', type = str, default = None, required=False,
                    help = 'dla catalogue')

parser.add_argument('--drq-qso-randoms', type = str, default = None, required=False, nargs='*',
                    help = 'qso randoms catalogues')

parser.add_argument('--drq-dla-randoms', type = str, default = None, required=False, nargs='*',
                    help = 'dla randoms catalogues')

parser.add_argument('--corr-dir', type = str, default = None, required=False,
                    help = 'directory to put correlations into')

parser.add_argument('--nproc', type = int, default = 64, required=False,
                    help = 'number of processes to use')

parser.add_argument('--nside', type = int, default = 16, required=False,
                    help = 'HEALPix nside for output files (must be 2^n)')

parser.add_argument('--run-all', action="store_true", default = False, required=False,
                    help = 'run all correlations')

parser.add_argument('--run-lyalya-lyalya', action="store_true", default = False, required=False,
                    help = 'run the auto correlation of lya absorption in lya region')

parser.add_argument('--run-lyalya-lyalyb', action="store_true", default = False, required=False,
                    help = 'run the cross correlation of lya absorption in lya and lyb regions')

parser.add_argument('--run-qso-qso', action="store_true", default = False, required=False,
                    help = 'run the auto correlation of qsos')

parser.add_argument('--run-dla-dla', action="store_true", default = False, required=False,
                    help = 'run the auto correlation of dlas')

parser.add_argument('--run-lyalya-qso', action="store_true", default = False, required=False,
                    help = 'run the cross correlation of lya absorption in lya region with qsos')

parser.add_argument('--run-lyalyb-qso', action="store_true", default = False, required=False,
                    help = 'run the cross correlation of lya absorption in lyb region with qsos')

parser.add_argument('--run-lyalya-dla', action="store_true", default = False, required=False,
                    help = 'run the cross correlation of lya absorption in lya region with dlas')

parser.add_argument('--run-lyalyb-dla', action="store_true", default = False, required=False,
                    help = 'run the cross correlation of lya absorption in lyb region with dlas')

parser.add_argument('--run-qso-dla', action="store_true", default = False, required=False,
                    help = 'run the qso-dla cross correlation')

parser.add_argument('--fid-Om', type=float, default=0.315, required=False,
                    help='Omega_matter(z=0) of fiducial LambdaCDM cosmology')

parser.add_argument('--fid-Or', type=float, default=0., required=False,
                    help='Omega_radiation(z=0) of fiducial LambdaCDM cosmology')

parser.add_argument('--rej', type=float, default=0.99, required=False,
                    help='proportional rejection to use for computing distortion matrices')

parser.add_argument('--no-project', action="store_true", default = False, required=False,
                    help = 'use the "no-project" option for Lya forests')

parser.add_argument('--no-dmats', action="store_true", default = False, required=False,
                    help = 'don\'t make dmats alongside correlations')

parser.add_argument('--no-remove-mean-lambda-obs', action="store_true", default = False, required=False,
                    help = 'use the "no-remove-mean-lambda-obs" option for Lya forest cross correlations')

parser.add_argument('--no-submit', action="store_true", default = False, required=False,
                    help = 'make the run scripts but do not submit the jobs')

args = parser.parse_args()

################################################################################

if args.run_all:
    args.run_lya_auto = True
    args.run_qso_auto = True
    args.run_dla_auto = True
    args.run_lya_aa_auto = True
    args.run_lya_qso_cross = True
    args.run_lya_dla_cross = True
    args.run_qso_dla_cross = True

################################################################################

#Functions to construct the job info dictionaries.
def make_cf_job_info(correlations_dir,corrname,zmin,zmax,deltas_dir,deltas_dir2=None,time=24,no_project=False,nodmat=False,rej=0.99):

    dir = os.path.join(correlations_dir,corrname)
    submit_utils.check_corr_dir(dir)
    out_dir = os.path.join(dir,'measurements')
    time = submit_utils.nh_to_hhmmss(time)

    header_info = {'queue':    'regular',
                   'time':     time,
                   'job_name': 'run_{}_{}_{}'.format(corrname,zmin,zmax),
                   'err_file': '{}_{}_{}_%j.err'.format(corrname,zmin,zmax),
                   'out_file': '{}_{}_{}_%j.out'.format(corrname,zmin,zmax),
                   }

    options = {'in-dir':        deltas_dir,
               'out':           os.path.join(out_dir,'cf_{}_{}_{}.fits.gz'.format(corrname,zmin,zmax)),
               'z-cut-min':     zmin,
               'z-cut-max':     zmax,
               }

    if deltas_dir2 is not None:
        options['in-dir2'] = deltas_dir2

    if no_project:
        options['no-project'] = ''

    job_info = {'dir':             dir,
                'header_info':     header_info,
                'picca_script':    'picca_cf.py',
                'options':         options,
                'run_script':      'run_{}_{}_{}.sl'.format(corrname,zmin,zmax),
                }

    if not nodmat:

        dmat_header_info = {'queue':    'regular',
                            'time':     time,
                            'job_name': 'run_dmat_{}_{}_{}'.format(corrname,zmin,zmax),
                            'err_file': 'dmat_{}_{}_{}_%j.err'.format(corrname,zmin,zmax),
                            'out_file': 'dmat_{}_{}_{}_%j.out'.format(corrname,zmin,zmax),
                            }

        dmat_options = {'in-dir':        deltas_dir,
                        'out':           os.path.join(out_dir,'dmat_{}_{}_{}.fits.gz'.format(corrname,zmin,zmax)),
                        'z-cut-min':     zmin,
                        'z-cut-max':     zmax,
                        'rej':           rej,
                        }

        if deltas_dir2 is not None:
            dmat_options['in-dir2'] = deltas_dir2

        if no_project:
            dmat_options['no-project'] = ''

        dmat_job_info = {'dir':             dir,
                         'header_info':     dmat_header_info,
                         'picca_script':    'picca_dmat.py',
                         'options':         dmat_options,
                         'run_script':      'run_dmat_{}_{}_{}.sl'.format(corrname,zmin,zmax),
                         }

    else:
        dmat_job_info = None

    return job_info, dmat_job_info

def make_co_job_info(correlations_dir,corrname,zmin,zmax,zcat,zcat2=None,time=24,zevolobj=1.44,zevolobj2=1.44):

    dir = os.path.join(correlations_dir,corrname)
    submit_utils.check_corr_dir(dir)
    out_dir = os.path.join(dir,'measurements')
    time = submit_utils.nh_to_hhmmss(time)

    header_info = {'queue':    'regular',
                   'time':     time,
                   'job_name': 'run_{}_{}_{}'.format(corrname,zmin,zmax),
                   'err_file': '{}_{}_{}_%j.err'.format(corrname,zmin,zmax),
                   'out_file': '{}_{}_{}_%j.out'.format(corrname,zmin,zmax),
                   }

    options = {'drq':           zcat,
               'z-evol-obj':    zevolobj,
               }

    if zcat2 is not None:
        options['drq2'] = zcat2
        options['z-evol-obj2'] = zevolobj2

    options = {**options,
               'out':           os.path.join(out_dir,'co_{}_{}_{}.fits.gz'.format(corrname,zmin,zmax)),
               'type-corr':     corr_type,
               'z-cut-min':     zmin,
               'z-cut-max':     zmax,
               }

    job_info = {'dir':             dir,
                'header_info':     header_info,
                'picca_script':    'picca_co.py',
                'options':         options,
                'run_script':      'run_{}_{}_{}.sl'.format(corrname,zmin,zmax),
                }

    return job_info

def make_xcf_job_info(correlations_dir,corrname,zmin,zmax,deltas_dir,zcat,time=24,no_project=False,no_remove_mean_lambda_obs=False,nodmat=False,rej=0.99,zevolobj=1.44):

    dir = os.path.join(correlations_dir,corrname)
    submit_utils.check_corr_dir(dir)
    out_dir = os.path.join(dir,'measurements')
    time = submit_utils.nh_to_hhmmss(time)

    header_info = {'queue':    'regular',
                   'time':     time,
                   'job_name': 'run_{}_{}_{}'.format(corrname,zmin,zmax),
                   'err_file': '{}_{}_{}_%j.err'.format(corrname,zmin,zmax),
                   'out_file': '{}_{}_{}_%j.out'.format(corrname,zmin,zmax),
                   }

    options = {'drq':           zcat,
               'in-dir':        deltas_dir,
               'out':           os.path.join(out_dir,'xcf_{}_{}_{}.fits.gz'.format(corrname,zmin,zmax)),
               'z-evol-obj':    zevolobj,
               'z-cut-min':     zmin,
               'z-cut-max':     zmax,
               }

    if no_project:
        options['no-project'] = ''
    if no_remove_mean_lambda_obs:
        options['no-remove-mean-lambda-obs'] = ''

    job_info = {'dir':          dir,
                'header_info':  header_info,
                'picca_script': 'picca_xcf.py',
                'options':      options,
                'run_script':   'run_{}_{}_{}.sl'.format(corrname,zmin,zmax),
                }

    if not nodmat:

        xdmat_header_info = {'queue':    'regular',
                             'time':     time,
                             'job_name': 'run_xdmat_{}_{}'.format(corrname,zmin,zmax),
                             'err_file': 'xdmat_{}_{}_{}_%j.err'.format(corrname,zmin,zmax),
                             'out_file': 'xdmat_{}_{}_{}_%j.out'.format(corrname,zmin,zmax),
                             }

        xdmat_options = {'drq':         zcat,
                         'in-dir':      deltas_dir,
                         'out':         os.path.join(out_dir,'xdmat_{}_{}_{}.fits.gz'.format(corrname,zmin,zmax)),
                         'z-evol-obj':  zevolobj,
                         'z-cut-min':   zmin,
                         'z-cut-max':   zmax,
                         'rej':         rej,
                         }

        if no_project:
            xdmat_options['no-project'] = ''

        xdmat_job_info = {'dir':             dir,
                          'header_info':     xdmat_header_info,
                          'picca_script':    'picca_xdmat.py',
                          'options':         xdmat_options,
                          'run_script':      'run_xdmat_{}_{}_{}.sl'.format(corrname,zmin,zmax),
                          }

    else:
        xdmat_job_info = None

    return job_info, xdmat_job_info


################################################################################

zbins = [(0.0,2.0),(2.0,2.2),(2.2,2.4),(2.4,2.6),(2.6,2.8),(2.8,3.0),(3.0,3.4),(3.4,10.0)]
global_job_info = {'options':   {'fid-Om':    args.fid_Om,
                                 'fid-Or':    args.fid_Or,
                                 'nside':     args.nside,
                                 'nproc':     args.nproc,
                                 },
                   }

job_time_dict = {'lyalya_lyalya':   {zbins[0]: 1.,   zbins[1]: 1.,  zbins[2]: 1.,  zbins[3]: 1.,
                                     zbins[4]: 1.,   zbins[5]: 1.,  zbins[6]: 1.,  zbins[7]: 1.,
                                     },
                 'lyalya_lyalyb':   {zbins[0]: 1.,   zbins[1]: 1.,  zbins[2]: 1.,  zbins[3]: 1.,
                                     zbins[4]: 1.,   zbins[5]: 1.,  zbins[6]: 1.,  zbins[7]: 1.,
                                     },
                 'qso_qso':         {'DD':   {zbins[0]: 3., zbins[1]: 3., zbins[2]: 3., zbins[3]: 3.,
                                              zbins[4]: 3., zbins[5]: 3., zbins[6]: 3., zbins[7]: 3.,
                                              },
                                     'RD':   {zbins[0]: 4., zbins[1]: 4., zbins[2]: 4., zbins[3]: 4.,
                                              zbins[4]: 4., zbins[5]: 4., zbins[6]: 4., zbins[7]: 4.,
                                              },
                                     'DR':   {zbins[0]: 4., zbins[1]: 4., zbins[2]: 4., zbins[3]: 4.,
                                              zbins[4]: 4., zbins[5]: 4., zbins[6]: 4., zbins[7]: 4.,
                                              },
                                     'RR':   {zbins[0]: 7., zbins[1]: 7., zbins[2]: 7., zbins[3]: 7.,
                                              zbins[4]: 7., zbins[5]: 7., zbins[6]: 7., zbins[7]: 7.,
                                              },
                                     },
                 'dla_dla':         {'DD':   {zbins[0]: 1., zbins[1]: 1., zbins[2]: 1., zbins[3]: 1.,
                                              zbins[4]: 1., zbins[5]: 1., zbins[6]: 1., zbins[7]: 1.,
                                              },
                                     'RD':   {zbins[0]: 2., zbins[1]: 2., zbins[2]: 2., zbins[3]: 2.,
                                              zbins[4]: 2., zbins[5]: 2., zbins[6]: 2., zbins[7]: 2.,
                                              },
                                     'DR':   {zbins[0]: 2., zbins[1]: 2., zbins[2]: 2., zbins[3]: 2.,
                                              zbins[4]: 2., zbins[5]: 2., zbins[6]: 2., zbins[7]: 2.,
                                              },
                                     'RR':   {zbins[0]: 4., zbins[1]: 4., zbins[2]: 4., zbins[3]: 4.,
                                              zbins[4]: 4., zbins[5]: 4., zbins[6]: 4., zbins[7]: 4.,
                                              },
                                     },
                 'lyalya_qso':      {'D':   {zbins[0]: 1., zbins[1]: 1., zbins[2]: 1., zbins[3]: 1.,
                                             zbins[4]: 1., zbins[5]: 1., zbins[6]: 1., zbins[7]: 1.,
                                             },
                                     'R':   {zbins[0]: 1., zbins[1]: 1., zbins[2]: 1., zbins[3]: 1.,
                                             zbins[4]: 1., zbins[5]: 1., zbins[6]: 1., zbins[7]: 1.,
                                             },
                                     },
                 'lyalya_dla':      {'D':   {zbins[0]: 1., zbins[1]: 1., zbins[2]: 1., zbins[3]: 1.,
                                             zbins[4]: 1., zbins[5]: 1., zbins[6]: 1., zbins[7]: 1.,
                                             },
                                     'R':   {zbins[0]: 1., zbins[1]: 1., zbins[2]: 1., zbins[3]: 1.,
                                             zbins[4]: 1., zbins[5]: 1., zbins[6]: 1., zbins[7]: 1.,
                                             },
                                     },
                 'lyalyb_qso':      {'D':   {zbins[0]: 1., zbins[1]: 1., zbins[2]: 1., zbins[3]: 1.,
                                             zbins[4]: 1., zbins[5]: 1., zbins[6]: 1., zbins[7]: 1.,
                                             },
                                     'R':   {zbins[0]: 1., zbins[1]: 1., zbins[2]: 1., zbins[3]: 1.,
                                             zbins[4]: 1., zbins[5]: 1., zbins[6]: 1., zbins[7]: 1.,
                                             },
                                     },
                 'lyalyb_dla':      {'D':   {zbins[0]: 1., zbins[1]: 1., zbins[2]: 1., zbins[3]: 1.,
                                             zbins[4]: 1., zbins[5]: 1., zbins[6]: 1., zbins[7]: 1.,
                                             },
                                     'R':   {zbins[0]: 1., zbins[1]: 1., zbins[2]: 1., zbins[3]: 1.,
                                             zbins[4]: 1., zbins[5]: 1., zbins[6]: 1., zbins[7]: 1.,
                                             },
                                     },
                 'qso_dla':   {'DD':   {zbins[0]: 2., zbins[1]: 2., zbins[2]: 2., zbins[3]: 2.,
                                              zbins[4]: 2., zbins[5]: 2., zbins[6]: 2., zbins[7]: 2.,
                                              },
                                     'RD':   {zbins[0]: 4., zbins[1]: 4., zbins[2]: 4., zbins[3]: 4.,
                                              zbins[4]: 4., zbins[5]: 4., zbins[6]: 4., zbins[7]: 4.,
                                              },
                                     'DR':   {zbins[0]: 4., zbins[1]: 4., zbins[2]: 4., zbins[3]: 4.,
                                              zbins[4]: 4., zbins[5]: 4., zbins[6]: 4., zbins[7]: 4.,
                                              },
                                     'RR':   {zbins[0]: 5., zbins[1]: 5., zbins[2]: 5., zbins[3]: 5.,
                                              zbins[4]: 5., zbins[5]: 5., zbins[6]: 5., zbins[7]: 5.,
                                              },
                                     },
                 }


################################################################################

njobs = 0
analysis_dir = submit_utils.AnalysisDir

for zbin in zbins:
    print(' -> looking at zbin {}'.format(zbin))
    zmin = zbin[0]
    zmax = zbin[1]

    ## Auto correlation of lya absorption in lya region
    if args.run_lyalya_lyalya:
        name = 'lyalya_lyalya'
        print(' -> -> looking at cf {}'.format(name))

        ## Make the job info
        time = job_time_dict[name][zbin]
        job_info, dmat_job_info = make_cf_job_info(args.corr_dir,name,zmin,zmax,args.deltas_dir_lya,time=time,no_project=args.no_project,nodmat=args.no_dmats,rej=args.rej)

        ## Submit the correlation job, plus the dmat if needed
        submit_utils.run_picca_job(job_info,global_job_info['options'],no_submit=args.no_submit)
        njobs += 1
        if not args.no_dmats:
            submit_utils.run_picca_job(dmat_job_info,global_job_info['options'],no_submit=args.no_submit)
            njobs += 1

    ## Cross correlation of lya absorption in lya and lyb regions
    if args.run_lyalya_lyalyb:
        name = 'lyalya_lyalyb'
        print(' -> -> looking at cf {}'.format(name))

        ## Make the job info
        time = job_time_dict[name][zbin]
        job_info, dmat_job_info = make_cf_job_info(args.corr_dir,name,zmin,zmax,args.deltas_dir_lya,deltas_dir2=args.deltas_dir_lyb,time=time,no_project=args.no_project,nodmat=args.no_dmats,rej=args.rej)

        ## Submit the correlation job, plus the dmat if needed
        submit_utils.run_picca_job(job_info,global_job_info['options'],no_submit=args.no_submit)
        njobs += 1
        if not args.no_dmats:
            submit_utils.run_picca_job(dmat_job_info,global_job_info['options'],no_submit=args.no_submit)
            njobs += 1

    ## Auto correlation of qsos
    if args.run_qso_qso:

        name = 'qso_qso'
        print(' -> -> looking at co {}'.format(name))

        ## Make the job info for DD
        tempname = name+'_DD'
        time = job_time_dict[name]['DD'][zbin]

        ## Submit the correlation job for DD
        job_info = make_co_job_info(args.corr_dir,tempname,zmin,zmax,args.drq_qso,time=time,zevolobj=1.44)
        submit_utils.run_picca_job(job_info,global_job_info['options'],no_submit=args.no_submit)
        njobs += 1

        if args.drq_qso_randoms is not None:

            print('WARN: randoms for the QSO auto are not yet fully implemented. Skipping...')

            """
            Really want to just do this from one catalogue?

            for i,randcat in enumerate(args.drq_qso_randoms):
                ## Make the job info for DD
                tempname = name+'_DR{}'.format(i)
                time = job_time_dict[name][zbin]

                ## Submit the correlation job for DD
                job_info = make_co_job_info(args.corr_dir,tempname,zmin,zmax,args.drq_qso,zcat2=randcat,time=time,zevolobj=1.44,zevolobj2=1.44)
                submit_utils.run_picca_job(job_info,global_job_info['options'],no_submit=args.no_submit)
                njobs += 1

                ## Make the job info for DD
                tempname = name+'_RR{}'.format(i)
                time = job_time_dict[name][zbin]

                ## Submit the correlation job for DD
                job_info = make_co_job_info(args.corr_dir,tempname,zmin,zmax,randcat,time=time,zevolobj=1.44)
                submit_utils.run_picca_job(job_info,global_job_info['options'],no_submit=args.no_submit)
                njobs += 1"""

    ## Auto correlation of dlas
    if args.run_dla_dla:

        name = 'dla_dla'
        print(' -> -> looking at co {}'.format(name))

        ## Make the job info for DD
        tempname = name+'_DD'
        time = job_time_dict[name]['DD'][zbin]

        ## Submit the correlation job for DD
        job_info = make_co_job_info(args.corr_dir,name,zmin,zmax,args.drq_dla,time=time,zevolobj=0.)
        submit_utils.run_picca_job(job_info,global_job_info['options'],no_submit=args.no_submit)
        njobs += 1

        if args.drq_dla_randoms is not None:

            print('WARN: randoms for the DLA auto are not yet fully implemented. Skipping...')

            """
            Really want to just do this from one catalogue?

            for i,randcat in enumerate(args.drq_dla_randoms):
                ## Make the job info for DD
                tempname = name+'_DR{}'.format(i)
                time = job_time_dict[name][zbin]

                ## Submit the correlation job for DD
                job_info = make_co_job_info(args.corr_dir,name,zmin,zmax,args.drq_dla,zcat2=randcat,time=time,zevolobj=0.,zevolobj2=0.)
                submit_utils.run_picca_job(job_info,global_job_info['options'],no_submit=args.no_submit)
                njobs += 1

                ## Make the job info for DD
                tempname = name+'_RR{}'.format(i)
                time = job_time_dict[name][zbin]

                ## Submit the correlation job for DD
                job_info = make_co_job_info(args.corr_dir,name,zmin,zmax,randcat,time=time,zevolobj=0.)
                submit_utils.run_picca_job(job_info,global_job_info['options'],no_submit=args.no_submit)
                njobs += 1"""

    ## Cross correlation of qsos and dlas
    if args.run_qso_dla:

        name = 'qso_dla'
        print(' -> -> looking at co {}'.format(name))

        ## Make the job info for DD
        tempname = name+'_DD'
        time = job_time_dict[name]['DD'][zbin]

        ## Submit the correlation job for DD
        job_info = make_co_job_info(args.corr_dir,name,zmin,zmax,args.drq_qso,zcat2=args.drq_dla,time=time,zevolobj=1.44,zevolobj2=0.)
        submit_utils.run_picca_job(job_info,global_job_info['options'],no_submit=args.no_submit)
        njobs += 1

        if args.drq_dla_randoms is not None:

            print('WARN: randoms for the QSO-DLA cross are not yet fully implemented. Skipping...')

            """
            WARN: This really isn't implemented at all, the below code is not correct

            for i,randcat in enumerate(args.drq_dla_randoms):
                ## Make the job info for DD
                tempname = name+'_DR{}'.format(i)
                time = job_time_dict[name][zbin]

                ## Submit the correlation job for DD
                job_info = make_co_job_info(args.corr_dir,name,zmin,zmax,args.drq_dla,zcat2=randcat,time=time,zevolobj=0.,zevolobj2=0.)
                submit_utils.run_picca_job(job_info,global_job_info['options'],no_submit=args.no_submit)
                njobs += 1

                ## Make the job info for DD
                tempname = name+'_RR{}'.format(i)
                time = job_time_dict[name][zbin]

                ## Submit the correlation job for DD
                job_info = make_co_job_info(args.corr_dir,name,zmin,zmax,randcat,time=time,zevolobj=0.)
                submit_utils.run_picca_job(job_info,global_job_info['options'],no_submit=args.no_submit)
                njobs += 1"""

    ## Cross correlation of lya absorption in lya region with qsos
    if args.run_lyalya_qso:

        name = 'lyalya_qso'
        print(' -> -> looking at xcf {}'.format(name))

        ## Make the job info for D
        tempname = name+'_D'
        time = job_time_dict[name]['D'][zbin]

        ## Submit the correlation job for D
        job_info, dmat_job_info = make_xcf_job_info(args.corr_dir,tempname,zmin,zmax,args.deltas_dir_lya,args.drq_qso,time=time,no_project=args.no_project,no_remove_mean_lambda_obs=args.no_remove_mean_lambda_obs,nodmat=args.no_dmats,rej=args.rej,zevolobj=1.44)
        submit_utils.run_picca_job(job_info,global_job_info['options'],no_submit=args.no_submit)
        njobs += 1
        if not args.no_dmats:
            submit_utils.run_picca_job(dmat_job_info,global_job_info['options'],no_submit=args.no_submit)
            njobs += 1

        if args.drq_qso_randoms is not None:

            for i,randcat in enumerate(args.drq_qso_randoms):
                ## Make the job info for R
                tempname = name+'_R{}'.format(i)
                time = job_time_dict[name]['R'][zbin]

                ## Submit the correlation job for R
                job_info, dmat_job_info = make_xcf_job_info(args.corr_dir,tempname,zmin,zmax,args.deltas_dir_lya,randcat,time=time,no_project=args.no_project,no_remove_mean_lambda_obs=args.no_remove_mean_lambda_obs,nodmat=args.no_dmats,rej=args.rej,zevolobj=1.44)
                submit_utils.run_picca_job(job_info,global_job_info['options'],no_submit=args.no_submit)
                njobs += 1
                if not args.no_dmats:
                    submit_utils.run_picca_job(dmat_job_info,global_job_info['options'],no_submit=args.no_submit)
                    njobs += 1

    ## Cross correlation of lya absorption in lyb region with qsos
    if args.run_lyalyb_qso:

        name = 'lyalyb_qso'
        print(' -> -> looking at xcf {}'.format(name))

        ## Make the job info for D
        tempname = name+'_D'
        time = job_time_dict[name]['D'][zbin]

        ## Submit the correlation job for D
        job_info, dmat_job_info = make_xcf_job_info(args.corr_dir,tempname,zmin,zmax,args.deltas_dir_lyb,args.drq_qso,time=time,no_project=args.no_project,no_remove_mean_lambda_obs=args.no_remove_mean_lambda_obs,nodmat=args.no_dmats,rej=args.rej,zevolobj=1.44)
        submit_utils.run_picca_job(job_info,global_job_info['options'],no_submit=args.no_submit)
        njobs += 1
        if not args.no_dmats:
            submit_utils.run_picca_job(dmat_job_info,global_job_info['options'],no_submit=args.no_submit)
            njobs += 1

        if args.drq_qso_randoms is not None:

            for i,randcat in enumerate(args.drq_qso_randoms):
                ## Make the job info for R
                tempname = name+'_R{}'.format(i)
                time = job_time_dict[name]['R'][zbin]

                ## Submit the correlation job for R
                job_info, dmat_job_info = make_xcf_job_info(args.corr_dir,tempname,zmin,zmax,args.deltas_dir_lyb,randcat,time=time,no_project=args.no_project,no_remove_mean_lambda_obs=args.no_remove_mean_lambda_obs,nodmat=args.no_dmats,rej=args.rej,zevolobj=1.44)
                submit_utils.run_picca_job(job_info,global_job_info['options'],no_submit=args.no_submit)
                njobs += 1
                if not args.no_dmats:
                    submit_utils.run_picca_job(dmat_job_info,global_job_info['options'],no_submit=args.no_submit)
                    njobs += 1

    ## Cross correlation of lya absorption in lya region with dlas
    if args.run_lyalya_dla:

        name = 'lyalya_dla'
        print(' -> -> looking at xcf {}'.format(name))

        ## Make the job info for D
        tempname = name+'_D'
        time = job_time_dict[name]['D'][zbin]

        ## Submit the correlation job for D
        job_info, dmat_job_info = make_xcf_job_info(args.corr_dir,tempname,zmin,zmax,args.deltas_dir_lya,args.drq_dla,time=time,no_project=args.no_project,no_remove_mean_lambda_obs=args.no_remove_mean_lambda_obs,nodmat=args.no_dmats,rej=args.rej,zevolobj=0.)
        submit_utils.run_picca_job(job_info,global_job_info['options'],no_submit=args.no_submit)
        njobs += 1
        if not args.no_dmats:
            submit_utils.run_picca_job(dmat_job_info,global_job_info['options'],no_submit=args.no_submit)
            njobs += 1

        if args.drq_dla_randoms is not None:

            for i,randcat in enumerate(args.drq_dla_randoms):
                ## Make the job info for R
                tempname = name+'_R{}'.format(i)
                time = job_time_dict[name]['R'][zbin]

                ## Submit the correlation job for R
                job_info, dmat_job_info = make_xcf_job_info(args.corr_dir,tempname,zmin,zmax,args.deltas_dir_lya,randcat,time=time,no_project=args.no_project,no_remove_mean_lambda_obs=args.no_remove_mean_lambda_obs,nodmat=args.no_dmats,rej=args.rej,zevolobj=0.)
                submit_utils.run_picca_job(job_info,global_job_info['options'],no_submit=args.no_submit)
                njobs += 1
                if not args.no_dmats:
                    submit_utils.run_picca_job(dmat_job_info,global_job_info['options'],no_submit=args.no_submit)
                    njobs += 1

    ## Cross correlation of lya absorption in lyb region with dlas
    if args.run_lyalyb_dla:

        name = 'lyalyb_dla'
        print(' -> -> looking at xcf {}'.format(name))

        ## Make the job info for D
        tempname = name+'_D'
        time = job_time_dict[name]['D'][zbin]

        ## Submit the correlation job for D
        job_info, dmat_job_info = make_xcf_job_info(args.corr_dir,tempname,zmin,zmax,args.deltas_dir_lyb,args.drq_dla,time=time,no_project=args.no_project,no_remove_mean_lambda_obs=args.no_remove_mean_lambda_obs,nodmat=args.no_dmats,rej=args.rej,zevolobj=0.)
        submit_utils.run_picca_job(job_info,global_job_info['options'],no_submit=args.no_submit)
        njobs += 1
        if not args.no_dmats:
            submit_utils.run_picca_job(dmat_job_info,global_job_info['options'],no_submit=args.no_submit)
            njobs += 1

        if args.drq_dla_randoms is not None:

            for i,randcat in enumerate(args.drq_dla_randoms):
                ## Make the job info for R
                tempname = name+'_R{}'.format(i)
                time = job_time_dict[name]['R'][zbin]

                ## Submit the correlation job for R
                job_info, dmat_job_info = make_xcf_job_info(args.corr_dir,tempname,zmin,zmax,args.deltas_dir_lyb,randcat,time=time,no_project=args.no_project,no_remove_mean_lambda_obs=args.no_remove_mean_lambda_obs,nodmat=args.no_dmats,rej=args.rej,zevolobj=0.)
                submit_utils.run_picca_job(job_info,global_job_info['options'],no_submit=args.no_submit)
                njobs += 1
                if not args.no_dmats:
                    submit_utils.run_picca_job(dmat_job_info,global_job_info['options'],no_submit=args.no_submit)
                    njobs += 1

    print('')


################################################################################

print('\nAll analyses for all realisations sent to the queue (total {} jobs).'.format(njobs))

################################################################################



"""#Functions to construct the job info dictionaries.
def make_lya_auto_job_info(meas_dir,zmin,zmax,deltas_dir,time=24,no_project=False,nodmat=False):

    dir = meas_dir + '/lya_auto/'
    out_dir = dir + '/correlations/'
    submit_utils.check_corr_dir(dir)
    time = submit_utils.nh_to_hhmmss(time)

    header_info = {'queue':    'regular',
                   'time':     time,
                   'job_name': 'run_lya_auto_{}_{}'.format(zmin,zmax),
                   'err_file': 'lya_auto_{}_{}_%j.err'.format(zmin,zmax),
                   'out_file': 'lya_auto_{}_{}_%j.out'.format(zmin,zmax),
                   }

    options = {'in-dir':        deltas_dir,
               'out':           out_dir+'/cf_lya_auto_{}_{}.fits.gz'.format(zmin,zmax),
               'z-cut-min':     zmin,
               'z-cut-max':     zmax,
               }

    if no_project:
        options['no-project'] = ''

    lya_auto_job_info = {'dir':             dir,
                         'header_info':     header_info,
                         'picca_script':    'picca_cf.py',
                         'options':         options,
                         'run_script':      'run_lya_auto_{}_{}.sh'.format(zmin,zmax),
                         }

    if not nodmat:

        dmat_header_info = {'queue':    'regular',
                            'time':     time,
                            'job_name': 'run_dmat_lya_auto_{}_{}'.format(zmin,zmax),
                            'err_file': 'dmat_lya_auto_{}_{}_%j.err'.format(zmin,zmax),
                            'out_file': 'dmat_lya_auto_{}_{}_%j.out'.format(zmin,zmax),
                            }

        dmat_options = {'in-dir':        deltas_dir,
                        'out':           out_dir+'/dmat_lya_auto_{}_{}.fits.gz'.format(zmin,zmax),
                        'z-cut-min':     zmin,
                        'z-cut-max':     zmax,
                        'rej':           0.99,
                        }

        if no_project:
            dmat_options['no-project'] = ''

        dmat_lya_auto_job_info = {'dir':             dir,
                                  'header_info':     dmat_header_info,
                                  'picca_script':    'picca_dmat.py',
                                  'options':         dmat_options,
                                  'run_script':      'run_dmat_lya_auto_{}_{}.sh'.format(zmin,zmax),
                                  }

    else:
        dmat_lya_auto_job_info = None

    return lya_auto_job_info, dmat_lya_auto_job_info

def make_qso_auto_job_info(meas_dir,zmin,zmax,zcat,zcat_rand,corr_type='DD',time=24):

    dir = meas_dir + '/qso_auto/'
    out_dir = dir + '/correlations/'
    submit_utils.check_corr_dir(dir)
    time = submit_utils.nh_to_hhmmss(time)

    header_info = {'queue':    'regular',
                   'time':     time,
                   'job_name': 'run_qso_auto_{}_{}_{}'.format(corr_type,zmin,zmax),
                   'err_file': 'qso_auto_{}_{}_{}_%j.err'.format(corr_type,zmin,zmax),
                   'out_file': 'qso_auto_{}_{}_{}_%j.out'.format(corr_type,zmin,zmax),
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

def make_dla_auto_job_info(meas_dir,zmin,zmax,zcat,zcat_rand,corr_type='DD',time=None):

    dir = meas_dir + '/dla_auto_lrmin1040.0_lrmax1200.0/'
    out_dir = dir + '/correlations/'
    submit_utils.check_corr_dir(dir)
    if time is None:
        time = 24
    time = submit_utils.nh_to_hhmmss(time)

    header_info = {'queue':    'regular',
                   'time':     time,
                   'job_name': 'run_dla_auto_{}_{}_{}'.format(corr_type,zmin,zmax),
                   'err_file': 'dla_auto_{}_{}_{}_%j.err'.format(corr_type,zmin,zmax),
                   'out_file': 'dla_auto_{}_{}_{}_%j.out'.format(corr_type,zmin,zmax),
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

def make_lya_qso_cross_job_info(meas_dir,zmin,zmax,deltas_dir,zcat,cat_type='D',time=24,no_project=False,no_remove_mean_lambda_obs=False,nodmat=False):

    dir = meas_dir + '/lya_qso_cross/'
    out_dir = dir + '/correlations/'
    submit_utils.check_corr_dir(dir)
    time = submit_utils.nh_to_hhmmss(time)

    header_info = {'queue':    'regular',
                   'time':     time,
                   'job_name': 'run_lya_qso_cross_{}_{}_{}'.format(cat_type,zmin,zmax),
                   'err_file': 'lya_qso_cross_{}_{}_{}_%j.err'.format(cat_type,zmin,zmax),
                   'out_file': 'lya_qso_cross_{}_{}_{}_%j.out'.format(cat_type,zmin,zmax),
                   }

    options = {'drq':                       zcat,
               'in-dir':                    deltas_dir,
               'out':                       out_dir+'xcf_lya_qso_cross_{}_{}_{}.fits.gz'.format(cat_type,zmin,zmax),
               'z-evol-obj':                1.44,
               'z-cut-min':                 zmin,
               'z-cut-max':                 zmax,
               }

    if no_project:
        options['no-project'] = ''
    if no_remove_mean_lambda_obs:
        options['no-remove-mean-lambda-obs'] = ''

    lya_qso_cross_job_info = {'dir':            dir,
                              'header_info':    header_info,
                              'picca_script':   'picca_xcf.py',
                              'options':        options,
                              'run_script':     'run_lya_qso_cross_{}_{}_{}.sh'.format(cat_type,zmin,zmax),
                              }

    if not nodmat:

        xdmat_header_info = {'queue':    'regular',
                             'time':     time,
                             'job_name': 'run_xdmat_lya_qso_cross_{}_{}'.format(zmin,zmax),
                             'err_file': 'xdmat_lya_qso_cross_{}_{}_%j.err'.format(zmin,zmax),
                             'out_file': 'xdmat_lya_qso_cross_{}_{}_%j.out'.format(zmin,zmax),
                             }

        xdmat_options = {'drq':                       zcat,
                         'in-dir':                    deltas_dir,
                         'out':                       out_dir+'xdmat_lya_qso_cross_{}_{}_{}.fits.gz'.format(cat_type,zmin,zmax),
                         'z-evol-obj':                1.44,
                         'z-cut-min':                 zmin,
                         'z-cut-max':                 zmax,
                         'rej':           0.99,
                         }

        if no_project:
            xdmat_options['no-project'] = ''

        xdmat_lya_qso_cross_job_info = {'dir':             dir,
                                        'header_info':     xdmat_header_info,
                                        'picca_script':    'picca_xdmat.py',
                                        'options':         xdmat_options,
                                        'run_script':      'run_xdmat_lya_qso_cross_{}_{}_{}.sh'.format(cat_type,zmin,zmax),
                                        }

    else:
        xdmat_lya_qso_cross_job_info = None

    return lya_qso_cross_job_info, xdmat_lya_qso_cross_job_info

def make_lya_dla_cross_job_info(meas_dir,zmin,zmax,deltas_dir,zcat,cat_type='D',time=24,no_project=False,no_remove_mean_lambda_obs=False,nodmat=False):

    dir = meas_dir + '/lya_dla_cross/'
    out_dir = dir + '/correlations/'
    submit_utils.check_corr_dir(dir)
    time = submit_utils.nh_to_hhmmss(time)

    header_info = {'queue':    'regular',
                   'time':     time,
                   'job_name': 'run_lya_dla_cross_{}_{}_{}'.format(cat_type,zmin,zmax),
                   'err_file': 'lya_dla_cross_{}_{}_{}_%j.err'.format(cat_type,zmin,zmax),
                   'out_file': 'lya_dla_cross_{}_{}_{}_%j.out'.format(cat_type,zmin,zmax),
                   }

    options = {'drq':                       zcat,
               'in-dir':                    deltas_dir,
               'out':                       out_dir+'xcf_lya_dla_cross_{}_{}_{}.fits.gz'.format(cat_type,zmin,zmax),
               'z-evol-obj':                0.0,
               'z-cut-min':                 zmin,
               'z-cut-max':                 zmax,
               }

    if no_project:
        options['no-project'] = ''
    if no_remove_mean_lambda_obs:
        options['no-remove-mean-lambda-obs'] = ''

    lya_dla_cross_job_info = {'dir':            dir,
                              'header_info':    header_info,
                              'picca_script':   'picca_xcf.py',
                              'options':        options,
                              'run_script':     'run_lya_dla_cross_{}_{}_{}.sh'.format(cat_type,zmin,zmax),
                              }

    if not nodmat:

        xdmat_header_info = {'queue':    'regular',
                             'time':     time,
                             'job_name': 'run_xdmat_lya_dla_cross_{}_{}'.format(zmin,zmax),
                             'err_file': 'xdmat_lya_dla_cross_{}_{}_%j.err'.format(zmin,zmax),
                             'out_file': 'xdmat_lya_dla_cross_{}_{}_%j.out'.format(zmin,zmax),
                             }

        xdmat_options = {'drq':                       zcat,
                         'in-dir':                    deltas_dir,
                         'out':                       out_dir+'xdmat_lya_dla_cross_{}_{}_{}.fits.gz'.format(cat_type,zmin,zmax),
                         'z-evol-obj':                0.0,
                         'z-cut-min':                 zmin,
                         'z-cut-max':                 zmax,
                         'rej':           0.99,
                         }

        if no_project:
            options['no-project'] = ''

        xdmat_lya_dla_cross_job_info = {'dir':             dir,
                                        'header_info':     xdmat_header_info,
                                        'picca_script':    'picca_xdmat.py',
                                        'options':         xdmat_options,
                                        'run_script':      'run_xdmat_lya_dla_cross_{}_{}_{}.sh'.format(cat_type,zmin,zmax),
                                        }

    else:
        xdmat_lya_dla_cross_job_info = None

    return lya_dla_cross_job_info, xdmat_lya_dla_cross_job_info

def make_qso_dla_cross_job_info(meas_dir,zmin,zmax,zcat_qso,zcat_qso_rand,zcat_dla,zcat_dla_rand,corr_type='DD',time=24):

    dir = meas_dir + '/qso_dla_cross/'
    out_dir = dir + '/correlations/'
    submit_utils.check_corr_dir(dir)
    time = submit_utils.nh_to_hhmmss(time)

    header_info = {'queue':    'regular',
                   'time':     time,
                   'job_name': 'run_qso_dla_cross_{}_{}_{}'.format(corr_type,zmin,zmax),
                   'err_file': 'qso_dla_cross_{}_{}_{}_%j.err'.format(corr_type,zmin,zmax),
                   'out_file': 'qso_dla_cross_{}_{}_{}_%j.out'.format(corr_type,zmin,zmax),
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







    if args.run_lya_auto:

        time = job_time_dict['lya_auto'][zbin]
        lya_auto_job_info, dmat_lya_auto_job_info = make_lya_auto_job_info(args.corr_dir,zmin,zmax,args.deltas_dir,time=time,no_project=args.no_project,nodmat=args.no_dmats)
        submit_utils.run_picca_job(lya_auto_job_info,global_job_info['options'],no_submit=args.no_submit)
        if not args.no_dmats:
            submit_utils.run_picca_job(dmat_lya_auto_job_info,global_job_info['options'],no_submit=args.no_submit)
        njobs += 1



    if args.run_lya_qso_cross or args.run_lya_dla_cross:
        print(' -> looking at zbin {}'.format(zbin))
        zmin = zbin[0]
        zmax = zbin[1]

    if args.run_lya_qso_cross:

        time = job_time_dict['lya_qso_cross']['D'][zbin]
        lya_qso_cross_job_info, xdmat_lya_qso_cross_job_info = make_lya_qso_cross_job_info(args.corr_dir,zmin,zmax,args.deltas_dir,args.drq_qso,cat_type='D',time=time,no_project=args.no_project,no_remove_mean_lambda_obs=args.no_remove_mean_lambda_obs,nodmat=args.no_dmats)
        submit_utils.run_picca_job(lya_qso_cross_job_info,global_job_info['options'],no_submit=args.no_submit)
        if not args.no_dmats:
            submit_utils.run_picca_job(xdmat_lya_qso_cross_job_info,global_job_info['options'],no_submit=args.no_submit)
        njobs += 1

        if args.drq_qso_randoms is not None:
            for i,randcat in enumerate(args.drq_qso_randoms):
                cat_type = 'R{}'.format(i)
                time = job_time_dict['lya_qso_cross']['R'][zbin]
                lya_qso_cross_job_info, xdmat_lya_qso_cross_job_info = make_lya_qso_cross_job_info(args.corr_dir,zmin,zmax,args.deltas_dir,randcat,cat_type=cat_type,time=time,no_project=args.no_project,no_remove_mean_lambda_obs=args.no_remove_mean_lambda_obs,nodmat=args.no_dmats)
                submit_utils.run_picca_job(lya_qso_cross_job_info,global_job_info['options'],no_submit=args.no_submit)
                if not args.no_dmats:
                    submit_utils.run_picca_job(xdmat_lya_qso_cross_job_info,global_job_info['options'],no_submit=args.no_submit)
                njobs += 1

    if args.run_lya_dla_cross:

        time = job_time_dict['lya_dla_cross']['D'][zbin]
        lya_dla_cross_job_info, xdmat_lya_dla_cross_job_info = make_lya_dla_cross_job_info(args.corr_dir,zmin,zmax,args.deltas_dir,args.drq_dla,cat_type='D',time=time,no_project=args.no_project,no_remove_mean_lambda_obs=args.no_remove_mean_lambda_obs,nodmat=args.no_dmats)
        submit_utils.run_picca_job(lya_dla_cross_job_info,global_job_info['options'],no_submit=args.no_submit)
        if not args.no_dmats:
            submit_utils.run_picca_job(xdmat_lya_dla_cross_job_info,global_job_info['options'],no_submit=args.no_submit)
        njobs += 1

        if args.drq_dla_randoms is not None:
            for i,randcat in enumerate(args.drq_dla_randoms):
                cat_type = 'R{}'.format(i)
                time = job_time_dict['lya_dla_cross']['R'][zbin]
                lya_dla_cross_job_info, xdmat_lya_dla_cross_job_info = make_lya_dla_cross_job_info(args.corr_dir,zmin,zmax,args.deltas_dir,randcat,cat_type=cat_type,time=time,no_project=args.no_project,no_remove_mean_lambda_obs=args.no_remove_mean_lambda_obs,nodmat=args.no_dmats)
                submit_utils.run_picca_job(lya_dla_cross_job_info,global_job_info['options'],no_submit=args.no_submit)
                if not args.no_dmats:
                    submit_utils.run_picca_job(xdmat_lya_dla_cross_job_info,global_job_info['options'],no_submit=args.no_submit)
                njobs += 1


    if args.run_qso_auto or args.run_dla_auto or args.run_qso_dla_cross:
        print(' -> looking at zbin {}'.format(zbin))
        zmin = zbin[0]
        zmax = zbin[1]

    if args.run_qso_auto:

        for corr_type in ['DD','RD','DR','RR']:
            time = job_time_dict['qso_auto'][corr_type][zbin]
            qso_auto_job_info = make_qso_auto_job_info(args.corr_dir,zmin,zmax,args.drq_qso,args.drq_qso_randoms,corr_type=corr_type,time=time)
            submit_utils.run_picca_job(qso_auto_job_info,global_job_info['options'],no_submit=args.no_submit)
            njobs += 1

    if args.run_dla_auto:

        for corr_type in ['DD','RD','DR','RR']:
            time = job_time_dict['dla_auto'][corr_type][zbin]
            dla_auto_job_info = make_dla_auto_job_info(args.corr_dir,zmin,zmax,args.drq_dla,args.drq_dla_randoms,corr_type=corr_type,time=time)
            submit_utils.run_picca_job(dla_auto_job_info,global_job_info['options'],no_submit=args.no_submit)
            njobs += 1

    if args.run_qso_dla_cross:

        for corr_type in ['DD','RD','DR','RR']:
            time = job_time_dict['qso_dla_cross'][corr_type][zbin]
            qso_dla_cross_job_info = make_qso_dla_cross_job_info(args.corr_dir,zmin,zmax,args.drq_qso,args.drq_qso_randoms,args.drq_dla,args.drq_dla_randoms,corr_type=corr_type,time=time)
            submit_utils.run_picca_job(qso_dla_cross_job_info,global_job_info['options'],no_submit=args.no_submit)
            njobs += 1
"""
