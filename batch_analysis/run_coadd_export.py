#!/usr/bin/env python

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

parser.add_argument('--v-realisations', type = int, default = [0], required=False,
                    help = 'realisation numbers of lyacolore realisations', nargs='*')

parser.add_argument('--nside', type = int, default = 16, required=False,
                    help = 'HEALPix nside for output files (must be 2^n)')

parser.add_argument('--export-individual-zbins', action="store_true", default = False, required=False,
                    help = 'export the qso-dla cross correlation')

parser.add_argument('--export-all', action="store_true", default = False, required=False,
                    help = 'export all correlations')

parser.add_argument('--export-lya-auto', action="store_true", default = False, required=False,
                    help = 'export the lya auto correlation')

parser.add_argument('--export-qso-auto', action="store_true", default = False, required=False,
                    help = 'export the qso auto correlation')

parser.add_argument('--export-dla-auto', action="store_true", default = False, required=False,
                    help = 'export the dla auto correlation')

parser.add_argument('--export-lya-aa-auto', action="store_true", default = False, required=False,
                    help = 'export the lya all absorber auto correlation')

parser.add_argument('--export-lya-qso-cross', action="store_true", default = False, required=False,
                    help = 'export the lya-qso cross correlation')

parser.add_argument('--export-lya-dla-cross', action="store_true", default = False, required=False,
                    help = 'export the lya-dla cross correlation')

parser.add_argument('--export-qso-dla-cross', action="store_true", default = False, required=False,
                    help = 'export the qso-dla cross correlation')

args = parser.parse_args()

################################################################################

a_dir = args.base_dir + '/analysis/'
if args.export_all:
    args.export_lya_auto = True
    args.export_qso_auto = True
    args.export_dla_auto = True
    args.export_lya_aa_auto = True
    args.export_lya_qso_cross = True
    args.export_lya_dla_cross = True
    args.export_qso_dla_cross = True

################################################################################

def coadd_export_cf(name,meas_dir,zbins):

    dir = meas_dir + '/' + name + '/'
    in_files = ''
    for zbin in zbins:
        in_file = dir + '/correlations/cf_{}_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
        if args.export_individual_zbins:
            out_file = dir + '/correlations/cf_exp_{}_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
            command = 'picca_export.py --data {} --out {}'.format(in_file,out_file)
            retcode = call(command,shell=True)
        in_files += in_file + ' '
    out_file = dir + '/correlations/cf_exp_{}.fits.gz'.format(name)
    command='picca_export_coadd_zint.py --data {} --out {} --no-dmat'.format(in_files,out_file)
    retcode = call(command,shell=True)

    return

def coadd_export_co(name,meas_dir,zbins):

    dir = meas_dir + '/' + name + '/'
    DD_files = ''
    DR_files = ''
    RD_files = ''
    RR_files = ''
    for zbin in zbins:
        DD_file = '/correlations/co_{}_DD_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
        DR_file = '/correlations/co_{}_DR_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
        RD_file = '/correlations/co_{}_RD_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
        RR_file = '/correlations/co_{}_RR_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
        if args.export_individual_zbins:
            out_file = dir + '/correlations/co_exp_{}_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
            command = 'picca_export_co.py --DD-file {} --DR-file {} --RD-file {} --RR-file {} --out {}'.format(DD_file,DR_file,RD_file,RR_file,out_file)
            retcode = call(command,shell=True)
        DD_files += DD_file + ' '
        DR_files += DR_file + ' '
        RD_files += RD_file + ' '
        RR_files += RR_file + ' '
    out_file = dir + '/correlations/co_exp_{}.fits.gz'.format(name)
    command='picca_export_coadd_zint_co.py --DD-files {} --DR-files {} --RD-files {} --RR-files {} --out {}'.format(DD_files,DR_files,RD_files,RR_files,out_file)
    retcode = call(command,shell=True)

    return

def coadd_export_xcf(name,meas_dir,zbins):

    dir = meas_dir + '/' + name + '/'
    D_files = ''
    R_files = ''
    for zbin in zbins:
        D_file = dir + '/correlations/xcf_{}_D_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
        R_file = dir + '/correlations/xcf_{}_R_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
        if args.export_individual_zbins:
            out_file = dir + '/correlations/xcf_exp_{}_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
            command = 'picca_export.py --data {} --out {} --remove-shuffled-correlation {}'.format(D_file,out_file,R_file)
            retcode = call(command,shell=True)
        D_files += D_file + ' '
        R_files += R_files + ' '
    out_file = dir + '/correlations/xcf_exp_{}.fits.gz'.format(name)
    command='picca_export_coadd_zint.py --data {} --out {} --no-dmat --remove-shuffled-correlation {}'.format(D_files,out_file,R_files)
    retcode = call(command,shell=True)

    return


"""
def coadd_export_lya_auto(meas_dir,zbins):

    dir = meas_dir + '/lya_auto/'
    in_files = ''
    for zbin in zbins:
        in_file = dir + '/correlations/cf_lya_auto_{}_{}.fits.gz'.format(zbin[0],zbin[1])
        if args.export_individual_zbins:
            out_file = dir + '/correlations/cf_exp_lya_auto_{}_{}.fits.gz'.format(zbin[0],zbin[1])
            command = 'picca_export.py --data {} --out {}'.format(in_file,out_file)
        in_files += in_file + ' '
    out_file = dir + '/correlations/cf_exp_lya_auto.fits.gz'
    command='picca_export_coadd_zint.py --data {} --out {} --no-dmat'.format(in_files,out_file)
    retcode = call(command,shell=True)

    return

def export_qso_auto(meas_dir,zbins):

    dir = meas_dir + '/qso_auto/'
    DD_files = ''
    DR_files = ''
    RD_files = ''
    RR_files = ''
    for zbin in zbins:
        DD_file =
        DD_files += dir + '/correlations/co_qso_auto_DD_{}_{}.fits.gz'.format(zbin[0],zbin[1]) + ' '
        DR_files += dir + '/correlations/co_qso_auto_DR_{}_{}.fits.gz'.format(zbin[0],zbin[1]) + ' '
        RD_files += dir + '/correlations/co_qso_auto_RD_{}_{}.fits.gz'.format(zbin[0],zbin[1]) + ' '
        RR_files += dir + '/correlations/co_qso_auto_RR_{}_{}.fits.gz'.format(zbin[0],zbin[1]) + ' '
    out_file = dir + '/correlations/co_exp_qso_auto.fits.gz'
    command='picca_export_coadd_zint_co.py --DD-files {} --DR-files {} --RD-files {} --RR-files {} --out {}'.format(DD_files,DR_files,RD_files,RR_files,out_file)
    retcode = call(command,shell=True)

    return

def export_dla_auto(meas_dir,zbins):

    dir = meas_dir + '/dla_auto/'
    DD_files = ''
    DR_files = ''
    RD_files = ''
    RR_files = ''
    for zbin in zbins:
        DD_files += dir + '/correlations/co_dla_auto_DD_{}_{}.fits.gz'.format(zbin[0],zbin[1]) + ' '
        DR_files += dir + '/correlations/co_dla_auto_DR_{}_{}.fits.gz'.format(zbin[0],zbin[1]) + ' '
        RD_files += dir + '/correlations/co_dla_auto_RD_{}_{}.fits.gz'.format(zbin[0],zbin[1]) + ' '
        RR_files += dir + '/correlations/co_dla_auto_RR_{}_{}.fits.gz'.format(zbin[0],zbin[1]) + ' '
    out_file = dir + '/correlations/co_exp_dla_auto.fits.gz'
    command='picca_export_coadd_zint_co.py --DD-files {} --DR-files {} --RD-files {} --RR-files {} --out {}'.format(DD_files,DR_files,RD_files,RR_files,out_file)
    retcode = call(command,shell=True)

    return

def coadd_export_lya_aa_auto(meas_dir,zbins):

    dir = meas_dir + '/lya_aa_auto/'
    in_files = ''
    for zbin in zbins:
        in_files += dir + '/correlations/cf_lya_aa_auto_{}_{}.fits.gz'.format(zbin[0],zbin[1]) + ' '
    out_file = dir + '/correlations/cf_exp_lya_aa_auto.fits.gz'
    command='picca_export_coadd_zint.py --data {} --out {} --no-dmat'.format(in_files,out_file)
    retcode = call(command,shell=True)

    return

def coadd_export_lya_qso_cross(meas_dir,zbins):

    dir = meas_dir + '/lya_qso_cross/'
    data_files = ''
    rand_files = ''
    for zbin in zbins:
        data_files += dir + '/correlations/xcf_lya_qso_cross_D_{}_{}.fits.gz'.format(zbin[0],zbin[1]) + ' '
        rand_files += dir + '/correlations/xcf_lya_qso_cross_R_{}_{}.fits.gz'.format(zbin[0],zbin[1]) + ' '
    out_file = dir + '/correlations/xcf_exp_lya_qso_cross.fits.gz'
    command='picca_export_coadd_zint.py --data {} --out {} --no-dmat --remove-shuffled-correlation {}'.format(data_files,out_file,rand_files)
    retcode = call(command,shell=True)

    return

def coadd_export_lya_dla_cross(meas_dir,zbins):

    dir = meas_dir + '/lya_dla_cross/'
    data_files = ''
    rand_files = ''
    for zbin in zbins:
        data_files += dir + '/correlations/xcf_lya_dla_cross_D_{}_{}.fits.gz'.format(zbin[0],zbin[1]) + ' '
        rand_files += dir + '/correlations/xcf_lya_dla_cross_R_{}_{}.fits.gz'.format(zbin[0],zbin[1]) + ' '
    out_file = dir + '/correlations/xcf_exp_lya_dla_cross.fits.gz'
    command='picca_export_coadd_zint.py --data {} --out {} --no-dmat --remove-shuffled-correlation {}'.format(data_files,out_file,rand_files)
    retcode = call(command,shell=True)

    return

def export_qso_dla_cross(meas_dir,zbins):

    dir = meas_dir + '/qso_dla_cross/'
    DD_files = ''
    DR_files = ''
    RD_files = ''
    RR_files = ''
    for zbin in zbins:
        DD_files += dir + '/correlations/co_qso_dla_cross_DD_{}_{}.fits.gz'.format(zbin[0],zbin[1]) + ' '
        DR_files += dir + '/correlations/co_qso_dla_cross_DR_{}_{}.fits.gz'.format(zbin[0],zbin[1]) + ' '
        RD_files += dir + '/correlations/co_qso_dla_cross_RD_{}_{}.fits.gz'.format(zbin[0],zbin[1]) + ' '
        RR_files += dir + '/correlations/co_qso_dla_cross_RR_{}_{}.fits.gz'.format(zbin[0],zbin[1]) + ' '
    out_file = dir + '/correlations/co_exp_qso_dla_cross.fits.gz'
    command='picca_export_coadd_zint_co.py --DD-files {} --DR-files {} --RD-files {} --RR-files {} --out {}'.format(DD_files,DR_files,RD_files,RR_files,out_file)
    retcode = call(command,shell=True)

    return
"""

################################################################################
#Merge the z bins together.

#cf_zbins = [(0.0,2.2),(2.2,2.6),(2.6,3.0),(3.0,10.0)]
cf_zbins = [(2.6,3.0),(3.0,10.0)]
xcf_zbins = cf_zbins
co_zbins = [(0.0,10.0)]

for v_rea in args.v_realisations:

    ver = 'v{}.{}.{}'.format(args.v_maj,args.v_min,v_rea)
    print('\nCoadding and exporting correlations for version {}:'.format(ver))

    #Check that the directories are constructed properly.
    ac_dir = a_dir+'/correlation_functions/'
    submit_utils.check_dir(ac_dir)
    acv_dir = ac_dir+'/'+ver+'/'
    submit_utils.check_dir(acv_dir)
    acvm_dir = acv_dir+'/measurements/'
    submit_utils.check_dir(acvm_dir)

    if args.export_lya_auto:
        coadd_export_cf('lya_auto',acvm_dir,cf_zbins)

    if args.export_qso_auto:
        coadd_export_co('qso_auto',acvm_dir,co_zbins)

    if args.export_dla_auto:
        coadd_export_co('dla_auto',acvm_dir,co_zbins)

    if args.export_lya_aa_auto:
        coadd_export_cf('lya_aa_auto',acvm_dir,cf_zbins)

    if args.export_lya_qso_cross:
        coadd_export_xcf('lya_qso_cross',acvm_dir,xcf_zbins)

    if args.export_lya_dla_cross:
        coadd_export_xcf('lya_dla_cross',acvm_dir,xcf_zbins)

    if args.export_qso_dla_cross:
        coadd_export_co('qso_dla_cross',acvm_dir,co_zbins)
