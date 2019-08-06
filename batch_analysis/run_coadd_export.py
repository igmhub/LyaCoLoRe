import numpy as np
from subprocess import call
import argparse
import os

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

def coadd_export_lya_auto(meas_dir,zbins):

    dir = meas_dir + '/lya_auto/'
    in_files = ''
    for zbin in zbins:
        in_files += dir + '/correlations/cf_lya_auto_{}_{}.fits.gz'.format(zbin[0],zbin[1]) + ' '
    out_file = dir + '/correlations/cf_exp_lya_auto.fits.gz'
    command='picca_export_coadd_zint.py --data {} --out {} --no-dmat'.format(in_files,out_file)
    retcode = call(command,shell=True)

    return

def export_qso_auto(meas_dir,zbins):

    dir = meas_dir + '/qso_auto/'

    for corr_type in ['DD','RD','DR','RR']:
        in_files = dir + '/correlations/co_qso_auto_{}_*_*.fits.gz'.format(corr_type)
        out_file = dir + '/correlations/co_qso_auto_{}.fits.gz'.format(corr_type)
        command='picca_export_coadd_zint.py --data {} --out {} --no-dmat'.format(in_files,out_file)
        retcode = call(command,shell=True)

    return

def export_dla_auto(meas_dir,zbins):

    dir = meas_dir + '/dla_auto/'

    for corr_type in ['DD','RD','DR','RR']:
        in_files = dir + '/correlations/co_dla_auto_{}_*_*.fits.gz'.format(corr_type)
        out_file = dir + '/correlations/co_dla_auto_{}.fits.gz'.format(corr_type)
        command='picca_export_coadd_zint.py --data {} --out {} --no-dmat'.format(in_files,out_file)
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
    for cat_type in ['D','R']:
        in_files = ''
        for zbin in zbins:
            in_files += dir + '/correlations/xcf_lya_qso_cross_{}_{}_{}.fits.gz'.format(cat_type,zbin[0],zbin[1]) + ' '
        out_file = dir + '/correlations/xcf_exp_lya_qso_cross.fits.gz'
    command='picca_export_coadd_zint.py --data {} --out {} --no-dmat'.format(in_files,out_file)
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
    command='picca_export_coadd_zint.py --data {} --out {} --no-dmat --remove-shuffled-correlation'.format(in_files,out_file)
    retcode = call(command,shell=True)

    return

def export_qso_dla_cross(meas_dir,zbins):

    dir = meas_dir + '/qso_dla_cross/'
    in_files = dir + '/correlations/cf_qso_dla_cross_*_*.fits.gz'
    out_file = dir + '/correlations/cf_qso_dla_cross.fits.gz'
    command='picca_export_coadd_zint.py --data {} --out {} --no-dmat'.format(in_files,out_file)
    retcode = call(command,shell=True)

    return

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
    #check_dir(ac_dir)
    acv_dir = ac_dir+'/'+ver+'/'
    #check_dir(acv_dir)
    acvm_dir = acv_dir+'/measurements/'
    #check_dir(acvm_dir)

    if args.export_lya_auto:
        coadd_export_lya_auto(acvm_dir,cf_zbins)

    if args.export_dla_auto:
        print('Cannot currently coadd co correlations: exporting each zbin separately.')
        export_dla_auto(acvm_dir,co_zbins)

    if args.export_qso_auto:
        print('Cannot currently coadd co correlations: exporting each zbin separately.')
        export_qso_auto(acvm_dir,co_zbins)

    if args.export_lya_aa_auto:
        coadd_export_lya_aa_auto(acvm_dir,cf_zbins)

    if args.export_lya_qso_cross:
        coadd_export_lya_qso_cross(acvm_dir,xcf_zbins)

    if args.export_lya_dla_cross:
        coadd_export_lya_dla_cross(acvm_dir,xcf_zbins)

    if args.export_qso_dla_cross:
        print('Cannot currently coadd co correlations: exporting each zbin separately.')
        export_qso_dla_cross(acvm_dir,co_zbins)

################################################################################
#Export everything.
"""
cf:
picca_export.py --data cf_0.5_200.fits.gz --out cf_exp_0.5_200.fits.gz

co:
picca_export_co.py --DD-file co_0.5_DD.fits.gz --RD-file co_0.5_RD.fits.gz --DR-file co_0.5_DR.fits.gz --RR-file co_0.5_RR.fits.gz --out co_exp_0.5.fits.gz

xcf:
picca_export.py --data xcf_0.5_200.fits.gz --out xcf_exp_0.5_200_shuffle.fits.gz --remove-shuffled-correlation xcf_0.5_200_shuffle.fits.gz
"""
