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
                    help = 'export the individual z bins')

parser.add_argument('--coadd-only', action="store_true", default = False, required=False,
                    help = 'coadd the (unexported) correlations')

parser.add_argument('--coadd-only-randoms', action="store_true", default = False, required=False,
                    help = 'coadd the (unexported) random correlations')

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

parser.add_argument('--stack-correlations', action="store_true", default = False, required=False,
                    help = 'stack the correlations')

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

    dir = meas_dir + '/' + name + '/correlations/'
    in_files = ''

    #Construct a list of the files which will be coadded.
    for zbin in zbins:
        in_file = dir + '/cf_{}_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
        in_files += in_file + ' '

        #If desired, export the zbins individually.
        if args.export_individual_zbins:
            out_file = dir + '/cf_exp_{}_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
            command = 'picca_export.py --data {} --out {}'.format(in_file,out_file)
            retcode = call(command,shell=True)

    #Construct the output names.
    out_file = dir + '/cf_exp_{}.fits.gz'.format(name)
    coadd_only_out_file = dir + '/cf_{}.fits.gz'.format(name)

    #Carry out the export/coadding.
    command='picca_export_coadd_zint.py --data {} --out {} --no-dmat --coadd-out {}'.format(in_files,out_file,coadd_only_out_file)
    retcode = call(command,shell=True)

    return

def coadd_export_xcf(name,meas_dir,zbins,coadd_only_randoms=False):

    dir = meas_dir + '/' + name + '/correlations/'
    D_files = ''
    R_files = ''
    for zbin in zbins:
        D_file = dir + '/xcf_{}_D_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
        R_file = dir + '/xcf_{}_R_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
        D_files += D_file + ' '
        R_files += R_file + ' '

        if args.export_individual_zbins:
            out_file = dir + '/xcf_exp_{}_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
            command = 'picca_export.py --data {} --out {} --remove-shuffled-correlation {}'.format(D_file,out_file,R_file)
            retcode = call(command,shell=True)

    #Construct the output names.
    out_file = dir + '/xcf_exp_{}.fits.gz'.format(name)
    coadd_only_out_file = dir + '/xcf_D_{}.fits.gz'.format(name)
    coadd_only_randoms_out_file = dir + '/xcf_R_{}.fits.gz'.format(name)

    #Carry out the export/coadding.
    command='picca_export_coadd_zint.py --data {} --out {} --no-dmat --remove-shuffled-correlation {} --coadd-out {} --coadd-out-shuffled {}'.format(D_files,out_file,R_files,coadd_only_out_file,coadd_only_randoms_out_file)
    retcode = call(command,shell=True)

    return

def coadd_export_co(name,meas_dir,zbins):

    dir = meas_dir + '/' + name + '/correlations/'
    DD_files = ''
    DR_files = ''
    RD_files = ''
    RR_files = ''

    #Construct lists of the files which will be coadded.
    for zbin in zbins:
        DD_file = dir + '/co_{}_DD_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
        DR_file = dir + '/co_{}_DR_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
        RD_file = dir + '/co_{}_RD_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
        RR_file = dir + '/co_{}_RR_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
        DD_files += DD_file + ' '
        DR_files += DR_file + ' '
        RD_files += RD_file + ' '
        RR_files += RR_file + ' '

        #If desired, export the zbins individually.
        if args.export_individual_zbins:
            out_file = dir + '/co_exp_{}_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
            command = 'picca_export_co.py --DD-file {} --DR-file {} --RD-file {} --RR-file {} --out {}'.format(DD_file,DR_file,RD_file,RR_file,out_file)
            retcode = call(command,shell=True)

    #Construct the output names.
    out_file = dir + '/co_exp_{}.fits.gz'.format(name)
    coadd_only_DD_out_file = dir + '/co_DD_{}.fits.gz'.format(name)
    coadd_only_DR_out_file = dir + '/co_DR_{}.fits.gz'.format(name)
    coadd_only_RD_out_file = dir + '/co_RD_{}.fits.gz'.format(name)
    coadd_only_RR_out_file = dir + '/co_RR_{}.fits.gz'.format(name)

    #Carry out the export/coadding.
    command='picca_export_coadd_zint_co.py --DD-files {} --DR-files {} --RD-files {} --RR-files {} --out {} --coadd-out-DD {} --coadd-out-DR {} --coadd-out-RD {} --coadd-out-RR {}'.format(DD_files,DR_files,RD_files,RR_files,out_file,coadd_only_DD_out_file,coadd_only_DR_out_file,coadd_only_RD_out_file,coadd_only_RR_out_file)
    retcode = call(command,shell=True)

    return

def stack_coadd_export_cf(name,corr_dir,vers,zbins):

    #Check that the relevant directories are set up correctly.
    s_dir = corr_dir + '/stack/'
    submit_utils.check_dir(s_dir)
    sm_dir = s_dir + '/measurements/'
    submit_utils.check_dir(sm_dir)
    smt_dir = sm_dir + name + '/'
    submit_utils.check_dir(smt_dir)
    smtc_dir = smt_dir + '/correlations/'
    submit_utils.check_dir(smtc_dir)

    #Stack the measurements by concatenating subsamples.
    fiin = []
    for ver in vers:
        fin = corr_dir + '/' + ver + '/measurements/' + name + '/correlations/cf_{}.fits.gz'.format(name)
        fiin += [fin]
    fout = smtc_dir + '/cf_{}.fits.gz'.format(name)
    submit_utils.concatenate_subsamples(fiin,fout,'cf')

    #Export the outfile.
    fout_exp = smtc_dir + '/cf_exp_{}.fits.gz'.format(name)
    command='picca_export.py --data {} --out {} --do-not-smooth-cov'.format(fout,fout_exp)
    retcode = call(command,shell=True)

    """
    in_files = ''
    for zbin in zbins:
        zbin_in_files = ''
        for in_dir in in_dirs:
            in_file = in_dir + '/cf_{}_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
            in_files += in_file + ' '
            zbin_in_files += in_file + ' '
        if args.export_individual_zbins:
            out_file = smtc_dir + '/cf_exp_{}_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
            command = 'picca_export_coadd_zint_co.py --data {} --out {} --no-dmat'.format(zbin_in_files,out_file)
            retcode = call(command,shell=True)
    out_file = smtc_dir + '/cf_exp_{}.fits.gz'.format(name)
    if coadd_only:
        coadd_only_out_file = smtc_dir + '/cf_{}.fits.gz'.format(name)
    else:
        coadd_only_out_file = None
    command='picca_export_coadd_zint.py --data {} --out {} --no-dmat --coadd-out {}'.format(in_files,out_file,coadd_only_out_file)
    retcode = call(command,shell=True)
    """

    return

def stack_coadd_export_xcf(name,corr_dir,vers,zbins):

    s_dir = corr_dir + '/stack/'
    submit_utils.check_dir(s_dir)
    sm_dir = s_dir + '/measurements/'
    submit_utils.check_dir(sm_dir)
    smt_dir = sm_dir + name + '/'
    submit_utils.check_dir(smt_dir)
    smtc_dir = smt_dir + '/correlations/'
    submit_utils.check_dir(smtc_dir)

    #For the data and the randoms.
    fouts = {}
    for dt in ['D','R']:

        #Stack the measurements by concatenating subsamples.
        fiin = []
        for ver in vers:
            fin = corr_dir + '/' + ver + '/measurements/' + name + '/correlations/xcf_{}_{}.fits.gz'.format(dt,name)
            fii += [fin]
        fout = smtc_dir + '/xcf_{}_{}.fits.gz'.format(dt,name)
        fouts[dt] = fout
        submit_utils.concatenate_subsamples(fiin,fout,'xcf')

    #Export the outfile.
    fout_exp = smtc_dir + '/xcf_exp_{}.fits.gz'.format(name)
    command='picca_export.py --data {} --out {} --do-not-smooth-cov'.format(fouts['D'],fout_exp)
    retcode = call(command,shell=True)

    """
    in_dirs = []
    for ver in vers:
        in_dir = corr_dir + '/' + ver + '/measurements/' + name + '/correlations/'
        in_dirs += [in_dir]

    D_files = ''
    R_files = ''
    for zbin in zbins:
        zbin_D_files = ''
        zbin_R_files = ''
        for in_dir in in_dirs:
            D_file = in_dir + '/xcf_{}_D_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
            R_file = in_dir + '/xcf_{}_R_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
            D_files += D_file + ' '
            R_files += R_file + ' '
            zbin_D_files += D_file + ' '
            zbin_R_files += R_file + ' '
        if args.export_individual_zbins:
            out_file = smtc_dir + '/xcf_exp_{}_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
            command = 'picca_export_coadd_zint_co.py --data {} --out {} --no-dmat --remove-shuffled-correlation {}'.format(zbin_D_files,out_file,zbin_R_files)
            retcode = call(command,shell=True)
    out_file = smtc_dir + '/xcf_exp_{}.fits.gz'.format(name)
    if coadd_only:
        coadd_only_out_file = smtc_dir + '/xcf_{}.fits.gz'.format(name)
    else:
        coadd_only_out_file = None
    if coadd_only:
        coadd_only_randoms_out_file = smtc_dir + '/xcf_R_{}.fits.gz'.format(name)
    else:
        coadd_only_randoms_out_file = None
    command='picca_export_coadd_zint.py --data {} --out {} --no-dmat --remove-shuffled-correlation {} --coadd-out {} --coadd-out-shuffled {}'.format(D_files,out_file,R_files,coadd_only_out_file,coadd_only_randoms_out_file)
    retcode = call(command,shell=True)
    """

    return

def stack_coadd_export_co(name,corr_dir,vers,zbins):


    s_dir = corr_dir + '/stack/'
    submit_utils.check_dir(s_dir)
    sm_dir = s_dir + '/measurements/'
    submit_utils.check_dir(sm_dir)
    smt_dir = sm_dir + name + '/'
    submit_utils.check_dir(smt_dir)
    smtc_dir = smt_dir + '/correlations/'
    submit_utils.check_dir(smtc_dir)

    #For the data and the randoms.
    fouts = {}
    for dt in ['DD','DR','RD','RR']:

        #Stack the measurements by concatenating subsamples.
        fiin = []
        for ver in vers:
            fin = corr_dir + '/' + ver + '/measurements/' + name + '/correlations/co_{}_{}.fits.gz'.format(dt,name)
            fii += [fin]
        fout = smtc_dir + '/co_{}_{}.fits.gz'.format(dt,name)
        fouts[dt] = fout
        submit_utils.concatenate_subsamples(fiin,fout,'xcf')

    #Export the outfile.
    fout_exp = smtc_dir + '/co_exp_{}.fits.gz'.format(name)
    command='picca_export_co.py --DD-file {} --DR-file {} --RD-file {} --RR-file {} --do-not-smooth-cov'.format(fouts['DD'],fouts['DR'],fouts['RD'],fouts['RR'],fout_exp)
    retcode = call(command,shell=True)

    """
    in_dirs = []
    for ver in vers:
        in_dir = corr_dir + '/' + ver + '/measurements/' + name + '/correlations/'
        in_dirs += [in_dir]

    DD_files = ''
    DR_files = ''
    RD_files = ''
    RR_files = ''
    for zbin in zbins:
        zbin_DD_files = ''
        zbin_DR_files = ''
        zbin_RD_files = ''
        zbin_RR_files = ''
        for in_dir in in_dirs:
            DD_file = in_dir + '/co_{}_DD_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
            DR_file = in_dir + '/co_{}_DR_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
            RD_file = in_dir + '/co_{}_RD_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
            RR_file = in_dir + '/co_{}_RR_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
            DD_files += DD_file + ' '
            DR_files += DR_file + ' '
            RD_files += RD_file + ' '
            RR_files += RR_file + ' '
            zbin_DD_files += DD_file + ' '
            zbin_DR_files += DR_file + ' '
            zbin_RD_files += RD_file + ' '
            zbin_RR_files += RR_file + ' '
        if args.export_individual_zbins:
            out_file = smtc_dir + '/co_exp_{}_{}_{}.fits.gz'.format(name,zbin[0],zbin[1])
            command = 'picca_export_coadd_zint_co.py --DD-files {} --DR-files {} --RD-files {} --RR-files {} --out {}'.format(zbin_DD_files,zbin_DR_files,zbin_RD_files,zbin_RR_files,out_file)
            retcode = call(command,shell=True)
    out_file = smtc_dir + '/co_exp_{}.fits.gz'.format(name)
    # TODO: implement the capability to coadd the unexported files.
    # Need to do so within picca first.
    command='picca_export_coadd_zint_co.py --DD-files {} --DR-files {} --RD-files {} --RR-files {} --out {}'.format(DD_files,DR_files,RD_files,RR_files,out_file)
    retcode = call(command,shell=True)
    """

    return


################################################################################
#Merge the z bins together.

cf_zbins = [(0.0,2.0),(2.0,2.2),(2.2,2.4),(2.4,2.6),(2.6,2.8),(2.8,3.0),(3.0,3.4),(3.4,10.0)]
xcf_zbins = cf_zbins
co_zbins = xcf_zbins

vers = []
for v_rea in args.v_realisations:

    ver = 'v{}.{}.{}'.format(args.v_maj,args.v_min,v_rea)
    vers += [ver]
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

if args.stack_correlations:
    stack_dir = ac_dir + '/stack/'
    submit_utils.check_corr_dir(stack_dir)

    if args.export_lya_auto:
        stack_coadd_export_cf('lya_auto',ac_dir,vers,cf_zbins)

    if args.export_qso_auto:
        stack_coadd_export_co('qso_auto',ac_dir,vers,cf_zbins)

    if args.export_dla_auto:
        stack_coadd_export_co('dla_auto',ac_dir,vers,cf_zbins)

    if args.export_lya_aa_auto:
        stack_coadd_export_cf('lya_aa_auto',ac_dir,vers,cf_zbins)

    if args.export_lya_qso_cross:
        stack_coadd_export_xcf('lya_qso_cross',ac_dir,vers,cf_zbins)

    if args.export_lya_dla_cross:
        stack_coadd_export_xcf('lya_dla_cross',ac_dir,vers,cf_zbins)

    if args.export_qso_dla_cross:
        stack_coadd_export_co('qso_dla_cross',ac_dir,vers,cf_zbins)
