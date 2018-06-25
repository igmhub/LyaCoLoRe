import numpy as np

import general

lya = 1215.67

#Function to extract RA values from a colore or picca format hdulist.
def get_RA(h,input_format):

    if input_format == 'physical_colore':
        RA = h[1].data['RA']
    elif input_format == 'gaussian_colore':
        RA = h[1].data['RA']
    elif input_format == 'picca':
        RA = h[3].data['RA']
    else:
        print('Error.')

    return RA

#Function to extract DEC values from a colore or picca format hdulist.
def get_DEC(h,input_format):

    if input_format == 'physical_colore':
        DEC = h[1].data['DEC']
    elif input_format == 'gaussian_colore':
        DEC = h[1].data['DEC']
    elif input_format == 'picca':
        DEC = h[3].data['DEC']
    else:
        print('Error.')

    return DEC

#Function to extract Z_QSO values from a colore or picca format hdulist.
def get_Z_QSO(h,input_format):

    if input_format == 'physical_colore':
        Z_QSO = h[1].data['Z_COSMO']
    elif input_format == 'gaussian_colore':
        Z_QSO = h[1].data['Z_COSMO']
    elif input_format == 'picca':
        Z_QSO = h[3].data['Z']
    else:
        print('Error.')

    return Z_QSO

#Function to extract DZ_RSD values from a colore.
def get_DZ_RSD(h,input_format):

    if input_format == 'physical_colore':
        DZ_RSD = h[1].data['DZ_RSD']
    elif input_format == 'gaussian_colore':
        DZ_RSD = h[1].data['DZ_RSD']
    elif input_format == 'picca':
        print('Error: DZ_RSD not stored in picca files.')
    else:
        print('Error.')

    return DZ_RSD

#Function to extract MOCKID values from a colore, picca or ID format hdulist.
def get_MOCKID(h,input_format,file_number):

    if input_format == 'physical_colore':
        #CoLoRe files do not have a MOCKID entry normally.
        #I am adding entries to any files processed via this code.
        #Hence we try to look for a MOCKID entry, and if this fails, we make one.
        try:
            MOCKID = h[1].data['MOCKID']
        except KeyError:
            h_N_qso = h[1].data.shape[0]
            row_numbers = list(range(h_N_qso))
            MOCKID = general.make_MOCKID(file_number,row_numbers)
    elif input_format == 'gaussian_colore':
        try:
            MOCKID = h[1].data['MOCKID']
        except KeyError:
            h_N_qso = h[1].data.shape[0]
            row_numbers = list(range(h_N_qso))
            MOCKID = general.make_MOCKID(file_number,row_numbers)
    elif input_format == 'picca':
        MOCKID = h[3].data['THING_ID']
    elif input_format == 'ID':
        MOCKID = h[1].data['MOCKID']

    return MOCKID

#Function to extract Z values from a colore or picca format hdulist.
def get_COSMO(h,input_format):

    lya = 1215.67

    if input_format == 'physical_colore':
        R = h[4].data['R']
        Z = h[4].data['Z']
        D = h[4].data['D']
        V = h[4].data['V']
    elif input_format == 'gaussian_colore':
        R = h[4].data['R']
        Z = h[4].data['Z']
        D = h[4].data['D']
        V = h[4].data['V']
    elif input_format == 'picca':
        LOGLAM_MAP = h[2].data
        Z = ((10**LOGLAM_MAP)/lya) - 1

        # TODO: Update this
        R = np.zeros(Z.shape)
        D = np.zeros(Z.shape)
        V = np.zeros(Z.shape)
    else:
        print('Error.')

    return R, Z, D, V

#Function to extract Z values from a colore or picca format hdulist.
def get_lya_lambdas(h,input_format):

    lya = 1215.67

    if input_format == 'physical_colore':
        Z = h[4].data['Z']
        lya_lambdas = lya*(1+Z)
    elif input_format == 'gaussian_colore':
        Z = h[4].data['Z']
        lya_lambdas = lya*(1+Z)
    elif input_format == 'picca':
        LOGLAM_MAP = h[2].data
        lya_lambdas = (10**LOGLAM_MAP)
    else:
        print('Error.')

    return lya_lambdas
