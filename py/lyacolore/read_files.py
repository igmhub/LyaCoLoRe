import numpy as np

from lyacolore import utils

lya = utils.lya_rest

## Function to get a stored value of SIGMA_G from an input file.
def get_SIGMA_G(h,file_format):

    if file_format == 'colore':
        try:
            SIGMA_G = h[4].header['SIGMA_G']
        except KeyError:
            SIGMA_G = None
    elif file_format == 'picca':
        print('INFO: picca does not store SIGMA_G.')
        SIGMA_G = None

    return SIGMA_G

## Function to get QSO information from an input file.
def get_QSO_data(h,file_format,wqso=None):

    TYPE = get_TYPE(h,file_format,wqso=wqso)
    RA = get_RA(h,file_format,wqso=wqso)
    DEC = get_DEC(h,file_format,wqso=wqso)
    Z_QSO = get_Z_QSO(h,file_format,wqso=wqso)
    DZ_RSD = get_DZ_RSD(h,file_format,wqso=wqso)

    return TYPE, RA, DEC, Z_QSO, DZ_RSD

#Function to extract TYPE values from a colore or picca format hdulist.
def get_TYPE(h,file_format,wqso=None):

    if file_format == 'colore':
        TYPE = h[1].data['TYPE']
    else:
        print('Error.')

    if wqso is not None:
        TYPE = TYPE[wqso]

    return TYPE

#Function to extract RA values from a colore or picca format hdulist.
def get_RA(h,file_format,wqso=None):

    if file_format == 'colore':
        RA = h[1].data['RA']
    elif file_format == 'picca':
        DEC = h[3].data['RA']
    else:
        print('Error.')

    if wqso is not None:
        RA = RA[wqso]

    return RA

#Function to extract DEC values from a colore or picca format hdulist.
def get_DEC(h,file_format,wqso=None):

    if file_format == 'colore':
        DEC = h[1].data['DEC']
    elif file_format == 'picca':
        DEC = h[3].data['DEC']
    else:
        print('Error.')

    if wqso is not None:
        DEC = DEC[wqso]

    return DEC

#Function to extract Z_QSO values from a colore or picca format hdulist.
def get_Z_QSO(h,file_format,wqso=None):

    if file_format == 'colore':
        Z_QSO = h[1].data['Z_COSMO']
    elif file_format == 'picca':
        Z_QSO = h[3].data['Z']
    else:
        print('Error.')

    if wqso is not None:
        Z_QSO = Z_QSO[wqso]

    return Z_QSO

#Function to extract DZ_RSD values from a colore.
def get_DZ_RSD(h,file_format,wqso=None):

    if file_format == 'colore':
        DZ_RSD = h[1].data['DZ_RSD']
    elif file_format == 'picca':
        print('Error: DZ_RSD not stored in picca files.')
    else:
        print('Error.')

    if wqso is not None:
        DZ_RSD = DZ_RSD[wqso]

    return DZ_RSD

#Function to extract MOCKID values from a colore, picca or ID format hdulist.
def get_MOCKID(h,file_format,file_number,wqso=None):

    if file_format == 'colore':
        #CoLoRe files do not have a MOCKID entry normally.
        #I am adding entries to any files processed via this code.
        #Hence we try to look for a MOCKID entry, and if this fails, we make one.
        try:
            MOCKID = h[1].data['MOCKID']
        except KeyError:
            h_N_qso = h[1].data.shape[0]
            row_numbers = list(range(h_N_qso))
            MOCKID = utils.make_MOCKID(file_number,row_numbers)
    elif file_format == 'picca':
        MOCKID = h[3].data['THING_ID']
    elif file_format == 'ID':
        MOCKID = h[1].data['MOCKID']

    if wqso is not None:
        MOCKID = MOCKID[wqso]

    return MOCKID

#Function to extract Z values from a colore or picca format hdulist.
def get_COSMO(h,file_format,wcell=None):

    if file_format == 'colore':
        R = h[4].data['R']
        Z = h[4].data['Z']
        D = h[4].data['D']
        V = h[4].data['V']
    elif file_format == 'picca':
        LOGLAM_MAP = h[2].data
        Z = ((10**LOGLAM_MAP)/lya) - 1

        # TODO: Update this
        R = np.zeros(Z.shape)
        D = np.zeros(Z.shape)
        V = np.zeros(Z.shape)
    else:
        print('Error.')

    if wcell is not None:
        R = R[wcell]
        Z = Z[wcell]
        D = D[wcell]
        V = V[wcell]

    return R, Z, D, V

#Function to extract Z values from a colore or picca format hdulist.
def get_lya_lambdas(h,file_format,wcell=None):

    if file_format == 'colore':
        Z = h[4].data['Z']
        lya_lambdas = lya*(1+Z)
    elif file_format == 'picca':
        LOGLAM_MAP = h[2].data
        lya_lambdas = (10**LOGLAM_MAP)
    else:
        print('Error.')

    if wcell is not None:
        lya_lambdas = lya_lambdas[wcell]

    return lya_lambdas

def get_skewer_rows(h,file_format,wqso=None,wcell=None):

    if file_format == 'colore':
        skewers = h[2].data
    elif file_format == 'picca':
        print('Error.')
    else:
        print('Error.')

    if wqso is not None:
        skewers = skewers[wqso,:]

    if wcell is not None:
        skewers = skewers[:,wcell]

    return skewers

def get_velocity_rows(h,file_format,wqso=None,wcell=None):

    if file_format == 'colore':
        vel = h[3].data
    elif file_format == 'picca':
        print('INFO: picca files do not store velocities.')
        vel = None
    else:
        print('Error.')

    if wqso is not None:
        vel = vel[wqso,:]

    if wcell is not None:
        vel = vel[:,wcell]

    return vel
