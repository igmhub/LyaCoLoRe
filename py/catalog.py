import numpy as np
from astropy.io import fits

from lyacolore import utils, read_files

#Function to extract data suitable for making ID files from a set of colore or picca format files.
def get_ID_data(filename,file_number,file_format,skewer_type,N_side,minimum_z=0.0,downsampling=1.0,QSO_filter=None,pixel_list=None,seed=0):

    ID_data = []
    cosmology = []
    N_pixels = 12*N_side**2

    #Open the file and extract the angular coordinate data.
    h = fits.open(filename)

    #Extract the component parts of the master file's data from h.
    _, RA, DEC, Z_QSO_NO_RSD, DZ_RSD = read_files.get_QSO_data(h,file_format)
    MOCKID = read_files.get_MOCKID(h,file_format,file_number)
    h_R, h_Z, h_D, h_V = read_files.get_COSMO(h,file_format)
    h.close()

    #Construct the remaining component parts of the master file's data.
    pixel_ID = utils.make_pixel_ID(N_side,RA,DEC)
    file_numbers = file_number * np.ones(RA.shape)

    #Calculate Z_QSO_RSD.
    Z_QSO_RSD = Z_QSO_NO_RSD + DZ_RSD

    #Join the pieces of the ID_data together.
    ID_data = list(zip(RA,DEC,Z_QSO_NO_RSD,Z_QSO_RSD,MOCKID,pixel_ID,file_numbers))

    #Sort the MOCKIDs and pixel_IDs into the right order: first by pixel number, and then by MOCKID.
    #Also filter out the objects with Z_QSO<minimum_z.
    dtype = [('RA', 'd'), ('DEC', 'd'), ('Z_QSO_NO_RSD', 'd'), ('Z_QSO_RSD', 'd'), ('MOCKID', int), ('PIXNUM', int), ('FILENUM', int)]
    ID = np.array(ID_data, dtype=dtype)
    ID = ID[ID['Z_QSO_NO_RSD']>minimum_z]
    ID_sort = np.sort(ID, order=['PIXNUM','MOCKID'])

    #Filter out QSOs that we don't want, either from a list of pixels or a Boolean function on (RA,DEC).
    if QSO_filter is not None:
        QSOs = QSO_filter(ID_sort['RA'],ID_sort['DEC'])
        ID_sort = ID_sort[QSOs]
    if pixel_list is not None:
        for pix in set(pixel_ID):
            if pix not in pixel_list:
                ID_sort = ID_sort[ID_sort['PIXNUM'] != pix]

    #Construct the cosmology array.
    cosmology_data = list(zip(h_R,h_Z,h_D,h_V))
    dtype = [('R', 'd'), ('Z', 'd'), ('D', 'd'), ('V', 'd')]
    cosmology = np.array(cosmology_data,dtype=dtype)

    #If we have any QSOs left, continue. Otherwise, return empties.
    N_qso = ID_sort.shape[0]
    if N_qso > 0:
        #Downsample if we want.
        if downsampling < 1.0:
            final_N_qso = round(N_qso*downsampling)
            gen = np.random.RandomState(seed=seed)
            random_QSOs = np.sort(gen.choice(N_qso,size=final_N_qso,replace=False))
            ID_sort = ID_sort[random_QSOs]

        #Make file-pixel map element and MOCKID lookup.
        pixel_ID_set = list(sorted(set([pixel for pixel in ID_sort['PIXNUM'] if pixel>=0])))
        file_pixel_map_element = np.zeros(N_pixels)
        MOCKID_lookup_element = {}
        for pixel in pixel_ID_set:
            file_pixel_map_element[pixel] = 1
            MOCKID_pixel_list = [ID_sort['MOCKID'][i] for i in range(len(ID_sort['PIXNUM'])) if ID_sort['PIXNUM'][i]==pixel]
            MOCKID_lookup_element = {**MOCKID_lookup_element,**{(file_number,pixel):MOCKID_pixel_list}}
    else:
        file_pixel_map_element = np.zeros(N_pixels)
        MOCKID_lookup_element = {}

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
def write_ID(filename,N_side,ID_data,cosmology_data=None,overwrite=False):

    #Make an appropriate header.
    header = fits.Header()
    header['NSIDE'] = N_side

    #Make a primary HDU.
    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)

    #Make the data into tables.
    hdu_ID = fits.BinTableHDU.from_columns(ID_data,header=header,name='CATALOG')
    if cosmology_data is not None:
        hdu_cosmology = fits.BinTableHDU.from_columns(cosmology_data,header=header,name='COSMO')
        hdulist = fits.HDUList([prihdu,hdu_ID,hdu_cosmology])
    else:
        hdulist = fits.HDUList([prihdu,hdu_ID])

    #Make the .fits file.
    hdulist.writeto(filename,overwrite=overwrite)
    hdulist.close()

    return

#Function to make the drq files needed for picca xcf functions.
def write_DRQ(filename,RSD_option,ID_data,N_side,overwrite=False):

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
    dtype = [('RA','f4'),('DEC','f4'),('Z','f4'),('THING_ID',int),('MJD','f4'),('FIBERID',int),('PLATE',int),('PIXNUM',int)]
    DRQ_data = np.array(list(zip(RA,DEC,Z,THING_ID,MJD,FID,PLATE,PIXNUM)),dtype=dtype)

    #Make an appropriate header.
    header = fits.Header()
    header['NSIDE'] = N_side

    #Create a new master file, with the same filename concatenated with '_picca_' and the RSD option chosen.
    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)
    hdu_DRQ = fits.BinTableHDU(DRQ_data,header=header)

    hdulist = fits.HDUList([prihdu,hdu_DRQ])
    hdulist.writeto(filename,overwrite=overwrite)
    hdulist.close()

    return

#Function to make a MOCKID lookup from master data.
def make_MOCKID_lookup(master_filepath,pixels):

    master = fits.open(master_filepath)
    master_data = master[1].data
    master.close()

    #Make a MOCKID lookup.
    master_data_pixel_set = set(master_data['PIXNUM'])
    pixels_set = set(pixels)
    pixel_list = list(sorted(master_data_pixel_set.intersection(pixels_set)))

    #For each pixel, find which QSOs are in it, and which files they come from.
    MOCKID_lookup = {}
    for pixel in pixel_list:

        pixel_indices = (master_data['PIXNUM']==pixel)
        pixel_MOCKIDs = master_data['MOCKID'][pixel_indices]
        pixel_file_number_list = list(sorted(set(master_data['FILENUM'][pixel_indices])))

        #For each file, find which relevant MOCKIDs are in it, and add them to the lookup.
        for file_number in pixel_file_number_list:

            pixel_file_indices = ((master_data['FILENUM'][pixel_indices])==file_number)
            pixel_file_MOCKIDs = pixel_MOCKIDs[pixel_file_indices]
            MOCKID_lookup = {**MOCKID_lookup,**{(file_number,pixel):list(pixel_file_MOCKIDs)}}

    return MOCKID_lookup
