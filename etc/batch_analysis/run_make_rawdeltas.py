import os
import sys
import glob
import numpy as np
import fitsio
import healpy

from picca import converters
from picca.utils import userprint

#lyacolore_path = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/test_lya_qso_cross/LyaCoLoRe_normal/'
#rawdeltas_path = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/test_lya_qso_cross/rawdeltas_normal/'

lyacolore_path = '/global/cfs/projectdirs/desi/mocks/lya_forest//london/v9.0/v9.0.1/'
rawdeltas_path = '/global/cfs/projectdirs/desi/mocks/lya_forest/picca/london/v9.0/v9.0.1/desi-raw'

nproc=60

# New function based on picca.converters.desi_convert_dla.
def desi_convert_qso(in_path, out_path, randoms=False):
    """Convert a catalog of QSO from a DESI mock format to the format used by picca
    Args:
        in_path: string
            Full path filename containing the QSO catalogue
        out_path: string
            Full path filename where the fits DLA catalogue will be written to
    """
    from_desi_key_to_picca_key = {
        'RA': 'RA',
        'DEC': 'DEC',
        'Z': 'Z_QSO_RSD',
        'THING_ID': 'MOCKID',
        'PLATE': 'MOCKID',
        'MJD': 'MOCKID',
        'FIBERID': 'MOCKID'
    }
    if randoms:
        from_desi_key_to_picca_key['Z'] = 'Z'
    # read catalogue
    cat = {}
    hdul = fitsio.FITS(in_path)
    for key, value in from_desi_key_to_picca_key.items():
        cat[key] = hdul['CATALOG'][value][:]
    hdul.close()
    userprint(("INFO: Found {} quasars").format( np.unique(cat['THING_ID']).size))
    # sort by THING_ID
    w = np.argsort(cat['THING_ID'])
    for key in cat:
        cat[key] = cat[key][w]

    for key in ['RA', 'DEC']:
        cat[key] = cat[key].astype('float64')

    # save results
    results = fitsio.FITS(out_path, 'rw', clobber=True)
    cols = list(cat.values())
    names = list(cat)
    results.write(cols, names=names, extname='CAT')
    results.close()

def desi_convert_dla(in_path, out_path, randoms=False):
    """Convert a catalog of DLA from a DESI format to the format used by picca
    Args:
        in_path: string
            Full path filename containing the ASCII DLA catalogue
        out_path: string
            Full path filename where the fits DLA catalogue will be written to
    """
    from_desi_key_to_picca_key = {
        'RA': 'RA',
        'DEC': 'DEC',
        'Z': 'Z_DLA_RSD',
        'ZQSO': 'Z_QSO_RSD',
        'NHI': 'N_HI_DLA',
        'THING_ID': 'MOCKID',
        'DLAID': 'DLAID',
        'PLATE': 'MOCKID',
        'MJD': 'MOCKID',
        'FIBERID': 'MOCKID',
    }
    if randoms:
        from_desi_key_to_picca_key['Z'] = 'Z_DLA'
        from_desi_key_to_picca_key.pop('NHI')
    # read catalogue
    cat = {}
    hdul = fitsio.FITS(in_path)
    for key, value in from_desi_key_to_picca_key.items():
        cat[key] = hdul['DLACAT'][value][:]
    hdul.close()
    userprint(("INFO: Found {} DLA from {} "
               "quasars").format(cat['Z'].size,
                                 np.unique(cat['THING_ID']).size))
    # sort by THING_ID
    w = np.argsort(cat['THING_ID'])
    for key in cat:
        cat[key] = cat[key][w]

    for key in ['RA', 'DEC']:
        cat[key] = cat[key].astype('float64')

    # save results
    results = fitsio.FITS(out_path, 'rw', clobber=True)
    cols = list(cat.values())
    names = list(cat)
    results.write(cols, names=names, extname='DLACAT')
    results.close()

# Make the QSO DRQs
desi_convert_qso(lyacolore_path+'/master.fits', rawdeltas_path+'/drq_qso.fits')
desi_convert_qso(lyacolore_path+'/master_randoms.fits', rawdeltas_path+'/drq_qso_randoms.fits',randoms=True)

# Make the DLA DRQs
desi_convert_dla(lyacolore_path+'/master_DLA.fits', rawdeltas_path+'/drq_dla.fits')
desi_convert_dla(lyacolore_path+'/master_DLA_randoms.fits', rawdeltas_path+'/drq_dla_randoms.fits',randoms=True)

# Make the deltas
converters.desi_convert_transmission_to_delta_files(rawdeltas_path+'/drq_qso.fits',
                                         rawdeltas_path+'/deltas/',
                                         in_dir=lyacolore_path,
                                         lambda_min=3600.,
                                         lambda_max=5500.,
                                         lambda_min_rest_frame=1040.,
                                         lambda_max_rest_frame=1200.,
                                         delta_log_lambda=3.e-4,
                                         max_num_spec=None,
                                         nproc=nproc)

