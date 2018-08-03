import numpy as np
from astropy.io import fits

import utils
import read_files

#Function to extract data suitable for making ID files from a set of colore or picca format files.
def get_ID_data(original_file_location,original_filename_structure,file_number,input_format,N_side,minimum_z=0.0):

    ID_data = []
    cosmology = []
    N_pixels = 12*N_side**2

    #Open the file and extract the angular coordinate data.
    filename = original_file_location + '/' + original_filename_structure.format(file_number)
    h = fits.open(filename)

    #Extract the component parts of the master file's data from h.
    RA = read_files.get_RA(h,input_format)
    DEC = read_files.get_DEC(h,input_format)
    Z_QSO_NO_RSD = read_files.get_Z_QSO(h,input_format)
    DZ_RSD = read_files.get_DZ_RSD(h,input_format)
    MOCKID = read_files.get_MOCKID(h,input_format,file_number)
    h_R, h_Z, h_D, h_V = read_files.get_COSMO(h,input_format)

    h.close()

    #Construct the remaining component parts of the master file's data.
    pixel_ID = utils.make_pixel_ID(N_side,RA,DEC)
    file_numbers = file_number * np.ones(RA.shape)

    #Calculate Z_QSO_RSD.
    Z_QSO_RSD = Z_QSO_NO_RSD + DZ_RSD

    #Join the pieces of the ID_data together.
    ID_data = list(zip(RA,DEC,Z_QSO_NO_RSD,Z_QSO_RSD,MOCKID,pixel_ID,file_numbers))

    #Sort the MOCKIDs and pixel_IDs into the right order: first by pixel number, and then by MOCKID.
    #Also filter out the objects with Z_QSO<minimum_z
    dtype = [('RA', 'd'), ('DEC', 'd'), ('Z_QSO_NO_RSD', 'd'), ('Z_QSO_RSD', 'd'), ('MOCKID', int), ('PIXNUM', int), ('FILENUM', int)]
    ID = np.array(ID_data, dtype=dtype)
    ID = ID[ID['Z_QSO_NO_RSD']>minimum_z]
    ID_sort = np.sort(ID, order=['PIXNUM','MOCKID'])

    #Make file-pixel map element and MOCKID lookup.
    pixel_ID_set = list(sorted(set([pixel for pixel in ID_sort['PIXNUM'] if pixel>=0])))
    file_pixel_map_element = np.zeros(N_pixels)
    MOCKID_lookup_element = {}
    for pixel in pixel_ID_set:
        file_pixel_map_element[pixel] = 1
        MOCKID_pixel_list = [ID_sort['MOCKID'][i] for i in range(len(ID_sort['PIXNUM'])) if ID_sort['PIXNUM'][i]==pixel]
        MOCKID_lookup_element = {**MOCKID_lookup_element,**{(file_number,pixel):MOCKID_pixel_list}}

    #Construct the cosmology array.
    cosmology_data = list(zip(h_R,h_Z,h_D,h_V))
    dtype = [('R', 'd'), ('Z', 'd'), ('D', 'd'), ('V', 'd')]
    cosmology = np.array(cosmology_data,dtype=dtype)

    return file_number, ID_sort, cosmology, file_pixel_map_element, MOCKID_lookup_element

#Function to join together the outputs from 'get_ID_data' in several multiprocessing processes.
def join_ID_data(results,N_side):

    file_numbers = []
    master_results = []
    bad_coordinates_results = []
    cosmology_results = []
    file_pixel_map_results = []
    MOCKID_lookup = {}

    for result in results:
        file_numbers += [result[0]]
        ID_result = result[1]
        master_results += [ID_result[ID_result['PIXNUM']>=0]]
        bad_coordinates_results += [ID_result[ID_result['PIXNUM']<0]]
        # TODO: Something to check that all cosmology results are the same
        cosmology_results = [result[2]]
        file_pixel_map_results += [result[3]]
        MOCKID_lookup = {**MOCKID_lookup,**result[4]}

    file_pixel_map = np.zeros((max(file_numbers)+1,12*(N_side**2)))
    for i, file_number in enumerate(file_numbers):
        file_pixel_map[file_number,:] = file_pixel_map_results[i]
    #print(master_results)
    master_data = np.concatenate(master_results)
    bad_coordinates_data = np.concatenate(bad_coordinates_results)
    cosmology_data = np.concatenate(cosmology_results)
    file_pixel_map = np.vstack(file_pixel_map_results)

    return master_data, bad_coordinates_data, cosmology_data, file_pixel_map, MOCKID_lookup

#Function to write a single ID file, given the data.
def write_ID(filename,ID_data,cosmology_data,N_side):

    #Make an appropriate header.
    header = fits.Header()
    header['NSIDE'] = N_side

    #Make the data into tables.
    hdu_ID = fits.BinTableHDU.from_columns(ID_data,header=header,name='CATALOG')
    hdu_cosmology = fits.BinTableHDU.from_columns(cosmology_data,header=header,name='COSMO')

    #Make a primary HDU.
    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)

    #Make the .fits file.
    hdulist = fits.HDUList([prihdu,hdu_ID,hdu_cosmology])
    hdulist.writeto(filename)
    hdulist.close()

    return

#Function to make the drq files needed for picca xcf functions.
def write_DRQ(filename,RSD_option,ID_data,N_side):

    #Extract data from the ID_data
    RA = ID_data['RA']
    DEC = ID_data['DEC']
    Z = ID_data['Z_QSO'+'_'+RSD_option]
    THING_ID = ID_data['MOCKID']
    PIXNUM = ID_data['PIXNUM']
    N_qso = RA.shape[0]

    #Fill in the data that is not specified.
    MJD = np.zeros(N_qso)
    FID = np.zeros(N_qso)
    PLATE = THING_ID

    #Make the data array.
    dtype = [('RA','f8'),('DEC','f8'),('Z','f8'),('THING_ID',int),('MJD','f8'),('FIBERID',int),('PLATE',int),('PIXNUM',int)]
    DRQ_data = np.array(list(zip(RA,DEC,Z,THING_ID,MJD,FID,PLATE,PIXNUM)),dtype=dtype)

    #Make an appropriate header.
    header = fits.Header()
    header['NSIDE'] = N_side

    #Create a new master file, with the same filename concatenated with '_picca_' and the RSD option chosen.
    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)
    hdu_DRQ = fits.BinTableHDU.from_columns(DRQ_data,header=header)

    hdulist = fits.HDUList([prihdu,hdu_DRQ])
    hdulist.writeto(filename)
    hdulist.close()

    return
