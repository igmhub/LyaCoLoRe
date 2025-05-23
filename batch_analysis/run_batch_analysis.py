#!/usr/bin/env python

import numpy as np
from subprocess import call
import argparse
import os

from lyacolore import submit_utils

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--base-dir', type = str, default = None, required=False,
                    help = 'directory containing the master file')

parser.add_argument('--v-maj', type = int, default = 9, required=False,
                    help = 'major version of lyacolore realisations')

parser.add_argument('--v-min', type = int, default = 0, required=False,
                    help = 'minor version of lyacolore realisations')

parser.add_argument('--v-realisations', type = int, default = 0, required=False,
                    help = 'realisation numbers of lyacolore realisations', nargs='*')

parser.add_argument('--nproc', type = int, default = 64, required=False,
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

#Functions to construct the job info dictionaries.
def make_lya_auto_job_info(meas_dir,ver,zmin,zmax,deltas_dir,time=None):

    dir = meas_dir + '/lya_auto/'
    out_dir = dir + '/correlations/'
    submit_utils.check_corr_dir(dir)
    if time is None:
        time = 24

    header_info = {'queue':    'regular',
                   'time':     '{}:00:00'.format(int(time)),
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

def make_qso_auto_job_info(meas_dir,ver,zmin,zmax,zcat,zcat_rand,corr_type='DD',time=None):

    dir = meas_dir + '/qso_auto/'
    out_dir = dir + '/correlations/'
    submit_utils.check_corr_dir(dir)
    if time is None:
        time = 24

    header_info = {'queue':    'regular',
                   'time':     '{}:00:00'.format(int(time)),
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

def make_dla_auto_job_info(meas_dir,ver,zmin,zmax,zcat,zcat_rand,corr_type='DD',time=None):

    dir = meas_dir + '/dla_auto_lrmin1040.0_lrmax1200.0/'
    out_dir = dir + '/correlations/'
    submit_utils.check_corr_dir(dir)
    if time is None:
        time = 24
    time = nh_to_hhmmss(time)

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

def make_lya_aa_auto_job_info(meas_dir,ver,zmin,zmax,deltas_dir,time=None):

    dir = meas_dir + '/lya_aa_auto/'
    out_dir = dir + '/correlations/'
    submit_utils.check_corr_dir(dir)
    if time is None:
        time = 24

    header_info = {'queue':    'regular',
                   'time':     '{}:00:00'.format(int(time)),
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

def make_lya_qso_cross_job_info(meas_dir,ver,zmin,zmax,deltas_dir,zcat,zcat_rand,cat_type='D',time=None):

    dir = meas_dir + '/lya_qso_cross/'
    out_dir = dir + '/correlations/'
    submit_utils.check_corr_dir(dir)
    if time is None:
        time = 24

    header_info = {'queue':    'regular',
                   'time':     '{}:00:00'.format(int(time)),
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

def make_lya_dla_cross_job_info(meas_dir,ver,zmin,zmax,deltas_dir,zcat,zcat_rand,cat_type='D',time=None):

    dir = meas_dir + '/lya_dla_cross/'
    out_dir = dir + '/correlations/'
    submit_utils.check_corr_dir(dir)
    if time is None:
        time = 24

    header_info = {'queue':    'regular',
                   'time':     '{}:00:00'.format(int(time)),
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

def make_qso_dla_cross_job_info(meas_dir,ver,zmin,zmax,zcat_qso,zcat_qso_rand,zcat_dla,zcat_dla_rand,corr_type='DD',time=None):

    dir = meas_dir + '/qso_dla_cross/'
    out_dir = dir + '/correlations/'
    submit_utils.check_corr_dir(dir)
    if time is None:
        time = 24

    header_info = {'queue':    'regular',
                   'time':     '{}:00:00'.format(int(time)),
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

################################################################################

cf_zbins = [(0.0,2.0),(2.0,2.2),(2.2,2.4),(2.4,2.6),(2.6,2.8),(2.8,3.0),(3.0,3.4),(3.4,10.0)]
xcf_zbins = cf_zbins
co_zbins = xcf_zbins
global_job_info = {'options':   {'fid-Om':    args.fid_Om,
                                 'fid-Or':    args.fid_Or,
                                 'nside':     args.nside,
                                 'nproc':     args.nproc,
                                 },
                   }

job_time_dict = {'lya_auto':        {cf_zbins[0]: 4.,   cf_zbins[1]: 5.,   cf_zbins[2]: 15.,  cf_zbins[3]: 15.,
                                     cf_zbins[4]: 10.,  cf_zbins[5]: 10.,  cf_zbins[6]: 3.,   cf_zbins[7]: 3.,
                                     },
                 'qso_auto':        {'DD':   {co_zbins[0]: 3., co_zbins[1]: 3., co_zbins[2]: 3., co_zbins[3]: 3.,
                                              co_zbins[4]: 3., co_zbins[5]: 3., co_zbins[6]: 3., co_zbins[7]: 3.,
                                              },
                                     'RD':   {co_zbins[0]: 6., co_zbins[1]: 6., co_zbins[2]: 6., co_zbins[3]: 6.,
                                              co_zbins[4]: 6., co_zbins[5]: 6., co_zbins[6]: 6., co_zbins[7]: 6.,
                                              },
                                     'DR':   {co_zbins[0]: 6., co_zbins[1]: 6., co_zbins[2]: 6., co_zbins[3]: 6.,
                                              co_zbins[4]: 6., co_zbins[5]: 6., co_zbins[6]: 6., co_zbins[7]: 6.,
                                              },
                                     'RR':   {co_zbins[0]: 7., co_zbins[1]: 7., co_zbins[2]: 7., co_zbins[3]: 7.,
                                              co_zbins[4]: 7., co_zbins[5]: 7., co_zbins[6]: 7., co_zbins[7]: 7.,
                                              },
                                     },
                 'dla_auto':        {'DD':   {co_zbins[0]: 1., co_zbins[1]: 1., co_zbins[2]: 1., co_zbins[3]: 1.,
                                              co_zbins[4]: 1., co_zbins[5]: 1., co_zbins[6]: 1., co_zbins[7]: 1.,
                                              },
                                     'RD':   {co_zbins[0]: 1., co_zbins[1]: 1., co_zbins[2]: 1., co_zbins[3]: 1.,
                                              co_zbins[4]: 1., co_zbins[5]: 1., co_zbins[6]: 1., co_zbins[7]: 1.,
                                              },
                                     'DR':   {co_zbins[0]: 1., co_zbins[1]: 1., co_zbins[2]: 1., co_zbins[3]: 1.,
                                              co_zbins[4]: 1., co_zbins[5]: 1., co_zbins[6]: 1., co_zbins[7]: 1.,
                                              },
                                     'RR':   {co_zbins[0]: 1.5, co_zbins[1]: 1.5, co_zbins[2]: 1.5, co_zbins[3]: 1.5,
                                              co_zbins[4]: 1.5, co_zbins[5]: 1.5, co_zbins[6]: 1.5, co_zbins[7]: 1.5,
                                              },
                                     },
                 'lya_aa_auto':     {cf_zbins[0]: 4.,   cf_zbins[1]: 5.,   cf_zbins[2]: 15.,  cf_zbins[3]: 15.,
                                     cf_zbins[4]: 10.,  cf_zbins[5]: 10.,  cf_zbins[6]: 3.,   cf_zbins[7]: 3.,
                                     },
                 'lya_qso_cross':   {'D':   {xcf_zbins[0]: 3., xcf_zbins[1]: 3., xcf_zbins[2]: 4., xcf_zbins[3]: 4.,
                                             xcf_zbins[4]: 4., xcf_zbins[5]: 4., xcf_zbins[6]: 3., xcf_zbins[7]: 3.,
                                             },
                                     'R':   {xcf_zbins[0]: 5., xcf_zbins[1]: 5., xcf_zbins[2]: 6., xcf_zbins[3]: 6.,
                                             xcf_zbins[4]: 6., xcf_zbins[5]: 6., xcf_zbins[6]: 5., xcf_zbins[7]: 5.,
                                             },
                                     },
                 'lya_dla_cross':   {'D':   {xcf_zbins[0]: 3., xcf_zbins[1]: 4., xcf_zbins[2]: 5., xcf_zbins[3]: 5.,
                                             xcf_zbins[4]: 4., xcf_zbins[5]: 3., xcf_zbins[6]: 3., xcf_zbins[7]: 3.,
                                             },
                                     'R':   {xcf_zbins[0]: 5., xcf_zbins[1]: 5., xcf_zbins[2]: 8., xcf_zbins[3]: 8.,
                                             xcf_zbins[4]: 5., xcf_zbins[5]: 5., xcf_zbins[6]: 3., xcf_zbins[7]: 3.,
                                             },
                                     },
                 'qso_dla_cross':   {'DD':   {co_zbins[0]: 2., co_zbins[1]: 2., co_zbins[2]: 2., co_zbins[3]: 2.,
                                              co_zbins[4]: 2., co_zbins[5]: 2., co_zbins[6]: 2., co_zbins[7]: 2.,
                                              },
                                     'RD':   {co_zbins[0]: 4., co_zbins[1]: 4., co_zbins[2]: 4., co_zbins[3]: 4.,
                                              co_zbins[4]: 4., co_zbins[5]: 4., co_zbins[6]: 4., co_zbins[7]: 4.,
                                              },
                                     'DR':   {co_zbins[0]: 4., co_zbins[1]: 4., co_zbins[2]: 4., co_zbins[3]: 4.,
                                              co_zbins[4]: 4., co_zbins[5]: 4., co_zbins[6]: 4., co_zbins[7]: 4.,
                                              },
                                     'RR':   {co_zbins[0]: 5., co_zbins[1]: 5., co_zbins[2]: 5., co_zbins[3]: 5.,
                                              co_zbins[4]: 5., co_zbins[5]: 5., co_zbins[6]: 5., co_zbins[7]: 5.,
                                              },
                                     },
                 }


################################################################################

njobs = 0
for v_rea in args.v_realisations:

    ver = 'v{}.{}.{}'.format(args.v_maj,args.v_min,v_rea)
    print('\nRunning analysis for version {}:'.format(ver))

    #Check that the directories are constructed properly.
    ac_dir = a_dir+'/correlation_functions/'
    submit_utils.check_dir(ac_dir)
    acv_dir = ac_dir+'/'+ver+'/'
    submit_utils.check_dir(acv_dir)
    acvm_dir = acv_dir+'/measurements/'
    submit_utils.check_dir(acvm_dir)

    #Define the location variables for this version.
    lya_deltas_loc = '{}/data/picca_input/{}/deltas_0.5/'.format(args.base_dir,ver)
    lya_aa_deltas_loc = '{}/data/picca_input/{}/deltas_0.5_Lyb_metals/'.format(args.base_dir,ver)
    zcat_qso_loc = '{}/data/picca_input/{}/zcat_0.5.fits'.format(args.base_dir,ver)
    zcat_dla_loc = '{}/data/picca_input/{}/zcat_DLA_0.5_lrmin1040.0_lrmax1200.0.fits'.format(args.base_dir,ver)
    zcat_qso_rand_loc = '{}/data/picca_input/{}/zcat_0.1_randoms.fits'.format(args.base_dir,ver)
    zcat_dla_rand_loc = '{}/data/picca_input/{}/zcat_DLA_0.1_randoms_lrmin1040.0_lrmax1200.0.fits'.format(args.base_dir,ver)

    for zbin in cf_zbins:

        print(' -> looking at zbin {}'.format(zbin))
        zmin = zbin[0]
        zmax = zbin[1]

        if args.run_lya_auto:

            time = job_time_dict['lya_auto'][zbin]
            lya_auto_job_info = make_lya_auto_job_info(acvm_dir,ver,zmin,zmax,lya_deltas_loc,time=time)
            submit_utils.run_picca_job(lya_auto_job_info,global_job_info['options'])
            njobs += 1

        if args.run_lya_aa_auto:

            time = job_time_dict['lya_aa_auto'][zbin]
            lya_aa_auto_job_info = make_lya_aa_auto_job_info(acvm_dir,ver,zmin,zmax,lya_aa_deltas_loc,time=time)
            submit_utils.run_picca_job(lya_aa_auto_job_info,global_job_info['options'])
            njobs += 1

    for zbin in xcf_zbins:

        print(' -> looking at zbin {}'.format(zbin))
        zmin = zbin[0]
        zmax = zbin[1]

        if args.run_lya_qso_cross:

            for cat_type in ['D','R']:
                time = job_time_dict['lya_qso_cross'][cat_type][zbin]
                lya_qso_cross_job_info = make_lya_qso_cross_job_info(acvm_dir,ver,zmin,zmax,lya_deltas_loc,zcat_qso_loc,zcat_qso_rand_loc,cat_type=cat_type,time=time)
                submit_utils.run_picca_job(lya_qso_cross_job_info,global_job_info['options'])
                njobs += 1

        if args.run_lya_dla_cross:

            for cat_type in ['D','R']:
                time = job_time_dict['lya_dla_cross'][cat_type][zbin]
                lya_dla_cross_job_info = make_lya_dla_cross_job_info(acvm_dir,ver,zmin,zmax,lya_deltas_loc,zcat_dla_loc,zcat_dla_rand_loc,cat_type=cat_type,time=time)
                submit_utils.run_picca_job(lya_dla_cross_job_info,global_job_info['options'])
                njobs += 1

    for zbin in co_zbins:

        print(' -> looking at zbin {}'.format(zbin))
        zmin = zbin[0]
        zmax = zbin[1]

        if args.run_qso_auto:

            for corr_type in ['DD','RD','DR','RR']:
                time = job_time_dict['qso_auto'][corr_type][zbin]
                qso_auto_job_info = make_qso_auto_job_info(acvm_dir,ver,zmin,zmax,zcat_qso_loc,zcat_qso_rand_loc,corr_type=corr_type,time=time)
                submit_utils.run_picca_job(qso_auto_job_info,global_job_info['options'])
                njobs += 1

        if args.run_dla_auto:

            for corr_type in ['DD','RD','DR','RR']:
                time = job_time_dict['dla_auto'][corr_type][zbin]
                dla_auto_job_info = make_dla_auto_job_info(acvm_dir,ver,zmin,zmax,zcat_dla_loc,zcat_dla_rand_loc,corr_type=corr_type,time=time)
                submit_utils.run_picca_job(dla_auto_job_info,global_job_info['options'])
                njobs += 1

        if args.run_qso_dla_cross:

            for corr_type in ['DD','RD','DR','RR']:
                time = job_time_dict['qso_dla_cross'][corr_type][zbin]
                qso_dla_cross_job_info = make_qso_dla_cross_job_info(acvm_dir,ver,zmin,zmax,zcat_qso_loc,zcat_qso_rand_loc,zcat_dla_loc,zcat_dla_rand_loc,corr_type=corr_type,time=time)
                submit_utils.run_picca_job(qso_dla_cross_job_info,global_job_info['options'])
                njobs += 1

################################################################################

print('\nAll analyses for all realisations sent to the queue (total {} jobs).'.format(njobs))

################################################################################
