import numpy as np
from astropy.io import fits
from astropy.table import Table
import healpy as hp
import os
import time

"""
To do:
 - add in options to trim pixels or make a random selection of quasars

 - inputs: input format, original_filename_structure, nodes, location, N_side, new_filename_structure, save_location, trim pixels?, reduce objects?, output format

 - sub-function possibilities

 - function to check whether the input format is correct

 - add bit to master to list files that went into it
"""

#Program to convert a .fits file of HEALPix pixels into several .fits files of smaller HEALPix pixels.
#Takes files of the format outputted by CoLoRe and converts into the same format that picca takes as input.

#Function to orchestrate the repixelisation process.
#First creates a master file recording the catalog of quasars being analysed.
#Then, creates per-HEALPix pixel files containing the density skewers in either colore or picca formats.
def split_pixels(original_location,original_filename_structure,file_numbers,input_format,N_side,save_location,new_filename_structure,output_format,z_min=0,existing_file_option=0):

    #Set constants.
    lya = 1215.67

    #Determine number of files and pixels.
    N_files = len(file_numbers)
    N_pixels = 12*N_side**2


    print('Creating catalog ID files.')

    #Get data for the master and 'bad coordinates' files, and a file-pixel map.
    #If file_pixel_map[i,j]=1, then file j contains objects in pixel i in it.
    master_data, bad_coordinates_data, file_pixel_map = get_ID_data(original_location,original_filename_structure,input_format,file_numbers,N_side,z_min)

    print('{} objects have been found with valid angular coordinates.'.format(master_data.shape[0]))
    print('{} objects have been found to have invalid angular coordinates.'.format(bad_coordinates_data.shape[0]))

    #Make a list of the pixels that the files cover.
    pixel_list = list(sorted(set(master_data['PIX'])))

    #Set up filenames.
    master_filename = save_location + '/' + 'zmin_{}_nside_{}_'.format(z_min,N_side) + 'master.fits'
    bad_coordinates_filename = save_location + '/' + 'zmin_{}_nside_{}_'.format(z_min,N_side) + 'bad_coordinates.fits'

    #Check for any pre-existing files with the filenames we want, and deal with them appropriately.
    #"Appropriately" is defined by 'existing_file_option', defined in 'repixelise.py'.
    master_filename, master_data = check_ID_filename(master_filename,master_data,existing_file_option)
    bad_coordinates_filename, bad_coordinates_data = check_ID_filename(bad_coordinates_filename,bad_coordinates_data,existing_file_option)

    #Write master file.
    write_ID(master_filename,master_data,N_side,z_min)
    print('Master file saved at :\n{}'.format(master_filename))

    #Write bad coordinates file.
    write_ID(bad_coordinates_filename,bad_coordinates_data,N_side,z_min)
    print('"Bad coordinates" file saved at:\n{}'.format(bad_coordinates_filename))

    print('Working on files named as:\n"{}",\nand labelled by:\nfile_numbers = {}\n'.format(original_filename_structure,file_numbers))
    print('Output files will be named as:\n"{}",\nand labelled by:\noutput format = {}\nzmin = {}\nnside = {}\npixels = {}\n'.format(new_filename_structure,output_format,z_min,N_side,pixel_list))

    #For each pixel, extract data, organise and save as a new file.
    for pixel in pixel_list:
        print('Working on HEALPix pixel number {} ({} of {})...'.format(pixel,pixel_list.index(pixel)+1,len(pixel_list)))

        #Determine which file numbers we need to look at for the current pixel.
        relevant_file_numbers = [file_number for file_number in file_numbers if file_pixel_map[pixel,file_number]==1]
        files_included = 0

        #For each relevant file, extract the data and aggregate over all files into a 'combined' object.
        for file_number in relevant_file_numbers:
            print(' -> Extracting data from file number {} ({} of {})...'.format(file_number,relevant_file_numbers.index(file_number)+1,len(relevant_file_numbers)))

            #Get the THING_IDs of the relevant quasars: those that are located in the current pixel, stored in the current file, and have z_qso<z_min.
            relevant_THING_IDs = [qso['THING_ID'] for qso in master_data if qso['PIX']==pixel and qso['FILE_NUMBER']==file_number]
            N_relevant_qso = len(relevant_THING_IDs)
            print('    -> {} relevant quasars found.'.format(N_relevant_qso))

            #If there are some relevant quasars, open the data file and make it into a simulation_data object.
            #We use simulation_data.get_reduced_data to avoid loading all of the file's data into the object.
            if N_relevant_qso > 0:
                filename = original_location + '/' + original_filename_structure.format(file_number)
                working = simulation_data.get_reduced_data(filename,file_number,input_format,relevant_THING_IDs,z_min)

            #Combine the data from the working file with that from the files already looked at.
            if files_included > 0:
                combined = simulation_data.combine_files(combined,working)
                files_included += 1
            else:
                combined = working
                files_included += 1

        print('Data extraction completed; {} quasars were found in total.'.format(combined.N_qso))

        #If we have extracted data fromn at least one file, we now try to save the 'combined' data as a new file.
        if files_included > 0:

            #Set up some useful headers.
            header = fits.Header()
            header['NSIDE'] = N_side
            header['PIX'] = pixel
            header['LYA'] = lya
            header['NQSO'] = combined.N_qso

            #Set up filename.
            filename = save_location + '/' + new_filename_structure.format(output_format,z_min,N_side,pixel)

            #Check for any issues with pre-existing files with the same filename.
            #Produce a final object and filename, as per 'existing_file_option'.
            final_filename, final_object = check_pixel_filename(filename,combined,output_format,pixel,existing_file_option,z_min)

            #Write the hdulist to the file format required.
            final_object.save(final_filename,header,output_format)
            print('Finished file written.\n')

    return

#Function to convert a list of numbers to a list of n-digit strings.
def numbers_to_strings(numbers,string_length):

    strings=[]
    for number in numbers:
        number_str = str(number)
        if len(number_str)<=string_length:
            string = '0'*(string_length-len(number_str))+number_str
        else:
            exit('The file number is too great to construct a unique THING_ID (more than 3 digits).')

        strings += [string]

    return strings

#Function to extract RA values from a colore or picca format hdulist.
def get_RA(h,input_format):

    if input_format == 'colore':
        RA = h[1].data['RA']
    elif input_format == 'picca':
        RA = h[3].data['RA']
    else:
        print('Error.')

    return RA

#Function to extract DEC values from a colore or picca format hdulist.
def get_DEC(h,input_format):

    if input_format == 'colore':
        DEC = h[1].data['DEC']
    elif input_format == 'picca':
        DEC = h[3].data['DEC']
    else:
        print('Error.')

    return DEC

#Function to extract Z_QSO values from a colore or picca format hdulist.
def get_Z_QSO(h,input_format):

    if input_format == 'colore':
        Z_QSO = h[1].data['Z_COSMO']
    elif input_format == 'picca':
        Z_QSO = h[3].data['Z']
    else:
        print('Error.')

    return Z_QSO

#Function to extract THING_ID values from a colore, picca or ID format hdulist.
def get_THING_ID(h,input_format,file_number):

    if input_format == 'colore':

        #CoLoRe files do not have a THING_ID entry normally.
        #I am adding entries to any files processed via this code.
        #Hence we try to look for a THING_ID entry, and if this fails, we make one.
        try:
            THING_ID = h[1].data['THING_ID']
        except KeyError:
            h_N_qso = h[1].data.shape[0]
            row_numbers = list(range(h_N_qso))
            THING_ID = make_THING_ID(file_number,row_numbers)

    elif input_format == 'picca':

        #Get THING_ID list.
        THING_ID = h[3].data['THING_ID']

    elif input_format == 'ID':

        THING_ID = h[1].data['THING_ID']

    return THING_ID

#Function to construct an array of THING_IDs given a file_number and a list of row_numbers.
def make_THING_ID(file_number,row_numbers):

    N_qso = len(row_numbers)
    node = '0'*(len(str(file_number))-3) + str(file_number)

    THING_ID = ['']*N_qso
    for i in range(N_qso):
        row_numbers[i] = str(row_numbers[i])
        if len(row_numbers[i])<=7:
            row_numbers[i] = '0'*(7-len(row_numbers[i]))+row_numbers[i]
        else:
            exit('The row number is too great to construct a unique THING_ID (more than 7 digits).')
        THING_ID[i] = int(node+row_numbers[i])

    THING_ID = np.array(THING_ID)

    return THING_ID

#Function to extract Z values from a colore or picca format hdulist.
def get_Z(h,input_format):

    lya = 1215.67

    if input_format == 'colore':
        Z = h[4].data['Z']
    elif input_format == 'picca':
        LOGLAM_MAP = h[2].data
        Z = ((10**LOGLAM_MAP)/lya) - 1
    else:
        print('Error.')

    return Z

#Function to determine in which HEALPix pixel each of a set of (RA,DEC) coordinates lies, given N_side.
def make_pixel_ID(N_side,RA,DEC):

    N_qso = RA.shape[0]

    #Convert DEC and RA in degrees to theta and phi in radians.
    theta = (np.pi/180.0)*(90.0-DEC)
    phi = (np.pi/180.0)*RA

    #Make a list of  the HEALPix pixel coordinate of each quasar.
    pixel_ID = ['']*N_qso
    for i in range(N_qso):
        #Check that the angular coordinates are valid. Put all objects with invalid coordinates into a non-realistic ID number (-1).
        if 0 <= theta[i] <= np.pi and 0 <= phi[i] <= 2*np.pi:
            pixel_ID[i] = int(hp.pixelfunc.ang2pix(N_side,theta[i],phi[i]))
        else:
            pixel_ID[i] = -1

    return pixel_ID

#Function to extract data suitable for making ID files from a set of colore or picca format files.
def get_ID_data(original_location,original_filename_structure,input_format,file_numbers,N_side,z_min):

    #Get list of nodes:
    #nodes_int = file_numbers

    #Convert the numbers into strings of length 3. e.g. 14 -> '014'
    #node_ID_length = 3
    #nodes = numbers_to_strings(nodes_int,node_ID_length)

    #Set up master_ID_data.
    #This will become a list of tuples, with one tuple for each qso.
    #Each tuple contains the qso's pixel, and its THING_ID.
    ID_data = []

    #Set up the file_pixel_map.
    #file_pixel_map[i,j] is 1 if the jth file contains qsos from the ith pixel, 0 otherwise.
    N_pixels = 12*N_side**2
    file_pixel_map = np.zeros((N_pixels,max(file_numbers)+1))

    for file_number in file_numbers:
        #Open the file and extract the angular coordinate data.
        filename = original_location + '/' + original_filename_structure.format(file_number)
        h = fits.open(filename)

        h_RA = get_RA(h,input_format)
        h_DEC = get_DEC(h,input_format)
        h_Z_QSO = get_Z_QSO(h,input_format)
        h_THING_ID = get_THING_ID(h,input_format,file_number)

        h_N_qso = h_RA.shape[0]

        h.close()

        if z_min > 0:
            indices = [i for i in range(h_N_qso) if h_Z_QSO[i]>z_min]
            RA = h_RA[indices]
            DEC = h_DEC[indices]
            THING_ID = h_THING_ID[indices]
        else:
            RA = h_RA
            DEC = h_DEC
            THING_ID = h_THING_ID

        pixel_ID = make_pixel_ID(N_side,RA,DEC)
        file_pixel_list = np.sort(list(set(pixel_ID)))
        N_qso = len(THING_ID)

        #Add the THING_ID and pixel_ID from this file to the ID list.
        for j in range(N_qso):
            ID_data += [(THING_ID[j],pixel_ID[j],file_number)]

        #Add information to file_pixel_map
        for pixel in file_pixel_list:
            if pixel >= 0:
                file_pixel_map[pixel,file_number]=1

    #Sort the THING_IDs and pixel_IDs into the right order: first by pixel number, and then by THING_ID.
    dtype = [('THING_ID', int), ('PIX', int), ('FILE_NUMBER', int)]
    ID = np.array(ID_data, dtype=dtype)
    ID_sort = np.sort(ID, order=['PIX','THING_ID'])

    #Separate the quasars with invalid coordinates from the ID data.
    ID_sort_filter = ID_sort[ID_sort['PIX']>=0]
    ID_sort_bad = ID_sort[ID_sort['PIX']<0]

    return ID_sort_filter, ID_sort_bad, file_pixel_map

#Function to write a single ID file, given the data.
def write_ID(filename,ID_data,N_side,z_min):

    #Add appropriate headers and make a table from the data.
    header = fits.Header()
    header['NSIDE'] = N_side
    bt = fits.BinTableHDU.from_columns(ID_data,header=header)

    #Make a primary HDU.
    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)

    #Make the .fits file.
    hdulist = fits.HDUList([prihdu,bt])
    hdulist.writeto(filename,overwrite=True)
    hdulist.close()

    return

#Function to check two lists or arrays for duplicate values.
def check_duplicates(A,B):

    sA = set(A)
    sB = set(B)
    if set.intersection(sA,sB) == set():
        result = False
    else:
        result = True

    return result

#Function to check whether an ID file with the proposed filename already exists and deal with it appropriately.
def check_ID_filename(proposed_filename,working_data,existing_file_option):

    #Set up check variables.
    working_filename = proposed_filename
    check = 0
    i = 2

    #Check to see if a file with this name already exists.
    #If so, deal with it in a sensible way. If not, move on.
    while check == 0 and os.path.exists(working_filename) == True:
        print('File named {} already exists.'.format(working_filename))

        #For this option, we always save the working data as a new file.
        if existing_file_option == 0:

            working_filename = proposed_filename[:-5] + '_v' + str(i) + '.fits'
            print('Will attempt to save new data in a new file named {}.'.format(working_filename))
            i = i+1

        #Check for duplicates: if there are duplicates, we will notice a change in the list lengths.
        else:
            #Produce lists of the THING_IDs in the 'working' dataset.
            working_THING_ID_list = list(working_data['THING_ID'])

            #Try to get a list of the THING_IDs in the 'existing' dataset.
            #We do not worry if the existing file has a different output_format.
            h = fits.open(working_filename)
            try:
                existing_THING_ID_list = list(get_THING_ID(h,'ID',''))
                format_match = 1
            except KeyError:
                print('Existing file is of a different format to requested "output_format".')
                format_match = 0
            h.close()

            duplicates_check = check_duplicates(existing_THING_ID_list,working_THING_ID_list)

            if format_match == 1:
                #If there are no duplicate THING_IDs.
                if duplicates_check == False:

                    #Attempt to merge it with the current file.
                    print('No duplicate entries found. Data will be added to existing file.')
                    h = fits.open(working_filename)
                    data_to_add = h[1].data
                    working_data = np.concatenate((working_data,data_to_add))
                    h.close()
                    check = 1

                #Else, if we chose option 1.
                elif existing_file_option == 1:

                    #Append the filename with '_i', and try to save again.
                    working_filename = proposed_filename[:-5] + '_v' + str(i) + '.fits'
                    print('Duplicate entries found. Will attempt to save new data in a new file named {}.'.format(working_filename))
                    i = i+1

                #Else, if we chose option 2.
                elif existing_file_option == 2:

                    print('Duplicate entries found. New data will be merged with existing data, omitting duplicates.')
                    #Merge the data with the current file, omitting duplicates.
                    THING_IDs_to_add_set = (set(working_THING_ID_list)-set(existing_THING_ID_list))
                    N_THING_IDs_to_add = len(THING_IDs_to_add_set)

                    print('{} objects will be added to existing file.'.format(N_THING_IDs_to_add))

                    rows = ['']*N_THING_IDs_to_add
                    j=0
                    for i, THING_ID in enumerate(working_THING_ID_list):
                        if THING_ID in THING_IDs_to_add_set:
                            rows[j] = i
                            j=j+1

                    data_to_add = working_data[rows]

                    h = fits.open(working_filename)
                    existing_data = h[1].data
                    h.close()

                    #data_to_add = [row for row in working_data if row['THING_ID'] in THING_IDs_to_add]
                    #data_to_add = np.array(data_to_add,dtype=existing_data.dtype)

                    working_data = np.concatenate((existing_data,data_to_add))
                    check = 1

            else:

                #Append the filename with '_i', and try to save again.
                working_filename = proposed_filename[:-5] + '_v' + str(i) + '.fits'
                print('Existing file is of a different format.')
                print('Will attempt to save new data in a new file named {}.'.format(proposed_filename))
                i = i+1

    final_filename = working_filename
    final_data = working_data

    return final_filename, final_data

#Function to check whether a per-pixel file with the proposed filename already exists and deal with it appropriately.
def check_pixel_filename(proposed_filename,working_object,output_format,pixel,existing_file_option,z_min):

    #Set up check variables.
    working_filename = proposed_filename
    check = 0
    i = 2

    #Check to see if a file with this name already exists.
    #If so, deal with it in a sensible way. If not, move on.
    while check == 0 and os.path.exists(working_filename) == True:
        print('File named {} already exists.'.format(working_filename))

        #For this option, we always save the working data as a new file.
        if existing_file_option == 0:

            working_filename = proposed_filename[:-5] + '_v' + str(i) + '.fits'
            print('Will attempt to save new data in a new file named {}.'.format(working_filename))
            i = i+1

        #Check for duplicates: if there are duplicates, we will notice a change in the list lengths.
        else:
            h = fits.open(working_filename)

            #Produce lists of the THING_IDs in the 'working' dataset.
            working_THING_ID_list = list(working_object.THING_ID)

            #Try to get a list of the THING_IDs in the 'existing' dataset.
            #We do not worry if the existing file has a different output_format.
            try:
                existing_THING_ID_list = list(get_THING_ID(h,output_format,pixel))
                format_match = 1
            except KeyError:
                print('Existing file is of a different format to requested "output_format".')
                format_match = 0

            h.close()

            duplicates_check = check_duplicates(existing_THING_ID_list,working_THING_ID_list)

            if format_match == 1:
                #If there are no duplicate THING_IDs.
                if duplicates_check == False:

                    #Attempt to merge it with the current file.
                    print('No duplicate entries found. Data will be added to existing file.')
                    existing = simulation_data.get_all_data(working_filename,pixel,output_format,z_min)
                    working_object = simulation_data.combine_files(existing,working_object)
                    check = 1

                #Else, if we chose option 1.
                elif existing_file_option == 1:

                    #Append the filename with '_i', and try to save again.
                    working_filename = proposed_filename[:-5] + '_v' + str(i) + '.fits'
                    print('Duplicate entries found. Will attempt to save new data in a new file named {}.'.format(working_filename))
                    i = i+1

                #Else, if we chose option 2.
                elif existing_file_option == 2:

                    print('Duplicate entries found. New data will be merged with existing data, omitting duplicates.')
                    #Merge the data with the current file, omitting duplicates.
                    THING_IDs_to_add_set = list(set(working_THING_ID_list)-set(existing_THING_ID_list))
                    N_THING_IDs_to_add = len(THING_IDs_to_add_set)

                    print('{} objects will be added to existing file.'.format(N_THING_IDs_to_add))

                    data_to_add = simulation_data.choose_qsos(working_object,THING_IDs_to_add_set)
                    existing_data = simulation_data.get_all_data(proposed_filename,pixel,output_format,z_min)

                    working_object = simulation_data.combine_files(existing_data,data_to_add)
                    check = 1

            else:

                #Append the filename with '_i', and try to save again.
                working_filename = proposed_filename[:-5] + '_v' + str(i) + '.fits'
                print('Existing file is of a different format.')
                print('Will attempt to save new data in a new file named {}.'.format(proposed_filename))
                i = i+1

    final_filename = working_filename
    final_object = working_object

    return final_filename, final_object

#Definition of the simulation_data class.
class simulation_data:
    #Initialisation function.
    def __init__(self,N_qso,N_cells,TYPE,RA,DEC,Z_QSO,DZ_RSD,THING_ID,PLATE,MJD,FIBER,DELTA_rows,VEL_rows,IVAR_rows,R,Z,D,V,LOGLAM_MAP):

        self.N_qso = N_qso
        self.N_cells = N_cells

        self.TYPE = TYPE
        self.RA = RA
        self.DEC = DEC
        self.Z_QSO = Z_QSO
        self.DZ_RSD = DZ_RSD
        self.THING_ID = THING_ID
        self.PLATE = PLATE
        self.MJD = MJD
        self.FIBER = FIBER

        self.DELTA_rows = DELTA_rows
        self.VEL_rows = VEL_rows
        self.IVAR_rows = IVAR_rows

        self.R = R
        self.Z = Z
        self.D = D
        self.V = V
        self.LOGLAM_MAP = LOGLAM_MAP

        return

    #Method to extract all data from an input file of a given format.
    @classmethod
    def get_all_data(cls,filename,file_number,input_format,z_min):

        lya = 1215.67
        h = fits.open(filename)

        h_Z = get_Z(h,input_format)

        #Calculate the first_relevant_cell.
        first_relevant_cell = np.argmax(h_Z >= z_min)

        if input_format == 'colore':

            #Extract data from the HDUlist.
            TYPE = h[1].data['TYPE']
            RA = h[1].data['RA']
            DEC = h[1].data['DEC']
            Z_QSO = h[1].data['Z_COSMO']
            DZ_RSD = h[1].data['DZ_RSD']
            DELTA_rows = h[2].data[:,first_relevant_cell:]
            VEL_rows = h[3].data[:,first_relevant_cell:]
            Z = h[4].data['Z'][first_relevant_cell:]
            R = h[4].data['R'][first_relevant_cell:]
            D = h[4].data['D'][first_relevant_cell:]
            V = h[4].data['V'][first_relevant_cell:]

            #Derive the number of quasars and cells in the file.
            N_qso = RA.shape[0]
            N_cells = Z.shape[0]

            #Derive the THING_ID and LOGLAM_MAP.
            THING_ID = get_THING_ID(h,input_format,file_number)
            LOGLAM_MAP = np.log10(lya*(1+Z))

            #Insert placeholder values for remaining variables.
            PLATE = np.zeros(N_qso)
            MJD = np.zeros(N_qso)
            FIBER = np.zeros(N_qso)
            IVAR_rows = np.ones((N_qso,N_cells))

        elif input_format == 'picca':

            #Extract data from the HDUlist.
            DELTA_rows = h[0].data.T[:,first_relevant_cell:]
            IVAR_rows = h[1].data.T[:,first_relevant_cell:]
            LOGLAM_MAP = h[2].data[first_relevant_cell:]
            RA = h[3].data['RA']
            DEC = h[3].data['DEC']
            Z_QSO = h[3].data['Z']
            PLATE = h[3].data['PLATE']
            MJD = h[3].data['MJD']
            FIBER = h[3].data['FIBER']
            THING_ID = h[3].data['THING_ID']

            #Derive the number of quasars and cells in the file.
            N_qso = RA.shape[0]
            N_cells = LOGLAM_MAP.shape[0]

            #Derive Z.
            Z = (10**LOGLAM_MAP)/lya - 1

            """
            Can we calculate DZ_RSD,R,D,V?
            """
            #Insert placeholder variables for remaining variables.
            TYPE = np.zeros(N_qso)
            DZ_RSD = np.zeros(N_qso)
            R = np.zeros(N_cells)
            D = np.zeros(N_cells)
            V = np.zeros(N_cells)
            VEL_rows = np.zeros((N_qso,N_cells))


        else:
            print('Input format not recognised: current options are "colore" and "picca".')
            print('Please choose one of these options and try again.')

        h.close()

        return cls(N_qso,N_cells,TYPE,RA,DEC,Z_QSO,DZ_RSD,THING_ID,PLATE,MJD,FIBER,DELTA_rows,VEL_rows,IVAR_rows,R,Z,D,V,LOGLAM_MAP)

    #NOT YET READY TO BE USED
    #Method to extract reduced data from a set of input files of a given format, with a given list of THING_IDs for each file.
    @classmethod
    def WIP(cls,file_infos,input_format,z_min):

        """
        NEED TO GENERATE FILE_INFOS TO PUT INTO THIS
        A DICTIONARY (?) OF FILENAME, FILE_NUMBER AND THING_IDS
        FORMATS STRING, INTEGER AND LIST
        """

        for file_info in file_infos:

            lya = 1215.67

            h = fits.open(file_info[filename])

            h_THING_ID = get_THING_ID(h,input_format,file_info[file_number])
            h_Z = get_Z(h,input_format)

            #Work out which rows in the hdulist we are interested in.
            rows = ['']*len(file_info[THING_IDs])
            s = set(file_info[THING_IDs])
            j = 0
            for i, qso in enumerate(h_THING_ID):
                if qso in s:
                    rows[j] = i
                    j = j+1

            #Calculate the first_relevant_cell.
            first_relevant_cell = np.argmax(h_Z >= z_min)

            TYPE = []
            RA = []
            DEC = []
            Z_QSO = []
            DZ_RSD = []
            DELTA_rows = []
            VEL_rows = []
            Z = []
            R = []
            D = []
            V = []
            N_qso = []
            N_cells = []
            THING_ID = []
            LOGLAM_MAP = []
            PLATE = []
            MJD = []
            FIBER = []
            IVAR_rows = []

            if input_format == 'colore':

                #Extract data from the HDUlist.
                TYPE = np.concatenate((TYPE,h[1].data['TYPE'][rows]))
                RA = np.concatenate((RA,h[1].data['RA'][rows]))
                DEC = np.concatenate((DEC,h[1].data['DEC'][rows]))
                Z_QSO = np.concatenate((Z_QSO,h[1].data['Z_COSMO'][rows]))
                DZ_RSD = np.concatenate((DZ_RSD,h[1].data['DZ_RSD'][rows]))

                DELTA_rows = np.concatenate((DELTA_rows,h[2].data[rows,first_relevant_cell:]),axis=0)

                VEL_rows = np.concatenate((VEL_rows,h[3].data[rows,first_relevant_cell:]),axis=0)

                Z = h[4].data['Z'][first_relevant_cell:]
                R = h[4].data['R'][first_relevant_cell:]
                D = h[4].data['D'][first_relevant_cell:]
                V = h[4].data['V'][first_relevant_cell:]

                #Derive the number of quasars and cells in the file.
                N_qso = N_qso + RA.shape[0]
                N_cells = Z.shape[0]

                #Derive the THING_ID and LOGLAM_MAP.
                THING_ID = np.concatenate((THING_ID,h_THING_ID[rows]))
                LOGLAM_MAP = np.log10(lya*(1+Z))

                #Insert placeholder values for remaining variables.
                PLATE = np.concatenate((PLATE,np.zeros(N_qso)))
                MJD = np.concatenate((MJD,np.zeros(N_qso)))
                FIBER = np.concatenate((FIBER,np.zeros(N_qso)))
                IVAR_rows = np.concatenate((IVAR_rows,np.ones((N_qso,N_cells))),axis=0)

            elif input_format == 'picca':

                """
                THIS NEEDS TO BE ADJUSTED IN THE SAME WAY THAT THE COLORE INPUT FORMAT HAS BEEN
                """

                #Extract data from the HDUlist.
                DELTA_rows = h[0].data.T[rows,first_relevant_cell:]

                IVAR_rows = h[1].data.T[rows,first_relevant_cell:]

                LOGLAM_MAP = h[2].data[first_relevant_cell:]

                RA = h[3].data['RA'][rows]
                DEC = h[3].data['DEC'][rows]
                Z_QSO = h[3].data['Z'][rows]
                PLATE = h[3].data['PLATE'][rows]
                MJD = h[3].data['MJD'][rows]
                FIBER = h[3].data['FIBER'][rows]
                THING_ID = h[3].data['THING_ID'][rows]

                #Derive the number of quasars and cells in the file.
                N_qso = RA.shape[0]
                N_cells = LOGLAM_MAP.shape[0]

                #Derive Z.
                Z = (10**LOGLAM_MAP)/lya - 1

                """
                Can we calculate DZ_RSD,R,D,V?
                """
                #Insert placeholder variables for remaining variables.
                TYPE = np.zeros(RA.shape[0])
                R = np.zeros(Z.shape[0])
                D = np.zeros(Z.shape[0])
                V = np.zeros(Z.shape[0])
                DZ_RSD = np.zeros(RA.shape[0])
                VEL_rows = np.zeros(DELTA_rows.shape)

            else:
                print('Input format not recognised: current options are "colore" and "picca".')
                print('Please choose one of these options and try again.')

            h.close()

        return cls(N_qso,N_cells,TYPE,RA,DEC,Z_QSO,DZ_RSD,THING_ID,PLATE,MJD,FIBER,DELTA_rows,VEL_rows,IVAR_rows,R,Z,D,V,LOGLAM_MAP)

    #Method to extract reduced data from an input file of a given format, with a given list of THING_IDs.
    @classmethod
    def get_reduced_data(cls,filename,file_number,input_format,THING_IDs,z_min):

        lya = 1215.67

        h = fits.open(filename)

        h_THING_ID = get_THING_ID(h,input_format,file_number)
        h_Z = get_Z(h,input_format)

        #Work out which rows in the hdulist we are interested in.
        rows = ['']*len(THING_IDs)
        s = set(THING_IDs)
        j = 0
        for i, qso in enumerate(h_THING_ID):
            if qso in s:
                rows[j] = i
                j = j+1

        #Calculate the first_relevant_cell.
        first_relevant_cell = np.argmax(h_Z >= z_min)

        if input_format == 'colore':

            #Extract data from the HDUlist.
            TYPE = h[1].data['TYPE'][rows]
            RA = h[1].data['RA'][rows]
            DEC = h[1].data['DEC'][rows]
            Z_QSO = h[1].data['Z_COSMO'][rows]
            DZ_RSD = h[1].data['DZ_RSD'][rows]

            DELTA_rows = h[2].data[rows,first_relevant_cell:]

            VEL_rows = h[3].data[rows,first_relevant_cell:]

            Z = h[4].data['Z'][first_relevant_cell:]
            R = h[4].data['R'][first_relevant_cell:]
            D = h[4].data['D'][first_relevant_cell:]
            V = h[4].data['V'][first_relevant_cell:]

            #Derive the number of quasars and cells in the file.
            N_qso = RA.shape[0]
            N_cells = Z.shape[0]

            #Derive the THING_ID and LOGLAM_MAP.
            THING_ID = h_THING_ID[rows]
            LOGLAM_MAP = np.log10(lya*(1+Z))

            #Insert placeholder values for remaining variables.
            PLATE = np.zeros(N_qso)
            MJD = np.zeros(N_qso)
            FIBER = np.zeros(N_qso)
            IVAR_rows = np.ones((N_qso,N_cells))

        elif input_format == 'picca':

            #Extract data from the HDUlist.
            DELTA_rows = h[0].data.T[rows,first_relevant_cell:]

            IVAR_rows = h[1].data.T[rows,first_relevant_cell:]

            LOGLAM_MAP = h[2].data[first_relevant_cell:]

            RA = h[3].data['RA'][rows]
            DEC = h[3].data['DEC'][rows]
            Z_QSO = h[3].data['Z'][rows]
            PLATE = h[3].data['PLATE'][rows]
            MJD = h[3].data['MJD'][rows]
            FIBER = h[3].data['FIBER'][rows]
            THING_ID = h[3].data['THING_ID'][rows]

            #Derive the number of quasars and cells in the file.
            N_qso = RA.shape[0]
            N_cells = LOGLAM_MAP.shape[0]

            #Derive Z.
            Z = (10**LOGLAM_MAP)/lya - 1

            """
            Can we calculate DZ_RSD,R,D,V?
            """
            #Insert placeholder variables for remaining variables.
            TYPE = np.zeros(RA.shape[0])
            R = np.zeros(Z.shape[0])
            D = np.zeros(Z.shape[0])
            V = np.zeros(Z.shape[0])
            DZ_RSD = np.zeros(RA.shape[0])
            VEL_rows = np.zeros(DELTA_rows.shape)

        else:
            print('Input format not recognised: current options are "colore" and "picca".')
            print('Please choose one of these options and try again.')

        h.close()

        return cls(N_qso,N_cells,TYPE,RA,DEC,Z_QSO,DZ_RSD,THING_ID,PLATE,MJD,FIBER,DELTA_rows,VEL_rows,IVAR_rows,R,Z,D,V,LOGLAM_MAP)

    @classmethod
    def combine_files(cls,A,B):

        N_qso = A.N_qso + B.N_qso
        """something to check N_cells is the same in both files"""
        N_cells = A.N_cells

        TYPE = np.concatenate((A.TYPE,B.TYPE),axis=0)
        RA = np.concatenate((A.RA,B.RA),axis=0)
        DEC = np.concatenate((A.DEC,B.DEC),axis=0)
        Z_QSO = np.concatenate((A.Z_QSO,B.Z_QSO),axis=0)
        DZ_RSD = np.concatenate((A.DZ_RSD,B.DZ_RSD),axis=0)
        THING_ID = np.concatenate((A.THING_ID,B.THING_ID),axis=0)
        PLATE = np.concatenate((A.PLATE,B.PLATE),axis=0)
        MJD = np.concatenate((A.MJD,B.MJD),axis=0)
        FIBER = np.concatenate((A.FIBER,B.FIBER),axis=0)

        DELTA_rows = np.concatenate((A.DELTA_rows,B.DELTA_rows),axis=0)
        VEL_rows = np.concatenate((A.VEL_rows,B.VEL_rows),axis=0)
        IVAR_rows = np.concatenate((A.IVAR_rows,B.IVAR_rows),axis=0)

        """
        Something to check this is ok?
        """
        Z = A.Z
        LOGLAM_MAP = A.LOGLAM_MAP
        R = A.R
        D = A.D
        V = A.V

        return cls(N_qso,N_cells,TYPE,RA,DEC,Z_QSO,DZ_RSD,THING_ID,PLATE,MJD,FIBER,DELTA_rows,VEL_rows,IVAR_rows,R,Z,D,V,LOGLAM_MAP)

    def save(self,filename,header,output_format):

        success = 0
        while success == 0:
            if output_format == 'colore':

                #Organise the data into colore-format arrays.
                colore_1_data = []
                for i in range(self.N_qso):
                    colore_1_data += [(self.TYPE[i],self.RA[i],self.DEC[i],self.Z_QSO[i],self.DZ_RSD[i],self.THING_ID[i])]

                dtype = [('TYPE', '>f4'), ('RA', '>f4'), ('DEC', '>f4'), ('Z_COSMO', '>f4'), ('DZ_RSD', '>f4'), ('THING_ID', int)]
                colore_1 = np.array(colore_1_data,dtype=dtype)

                colore_2 = self.DELTA_rows
                colore_3 = self.VEL_rows

                colore_4_data = []
                for i in range(self.N_cells):
                    colore_4_data += [(self.R[i],self.Z[i],self.D[i],self.V[i])]

                dtype = [('R', '>f4'), ('Z', '>f4'), ('D', '>f4'), ('V', '>f4')]
                colore_4 = np.array(colore_4_data,dtype=dtype)

                #Construct HDUs from the data arrays.
                prihdr = fits.Header()
                prihdu = fits.PrimaryHDU(header=prihdr)
                cols_CATALOG = fits.ColDefs(colore_1)
                hdu_CATALOG = fits.BinTableHDU.from_columns(cols_CATALOG,header=header,name='CATALOG')
                hdu_DELTA = fits.ImageHDU(data=colore_2,header=header,name='DELTA')
                hdu_VEL = fits.ImageHDU(data=colore_3,header=header,name='VELOCITY')
                cols_COSMO = fits.ColDefs(colore_4)
                hdu_COSMO = fits.BinTableHDU.from_columns(cols_COSMO,header=header,name='CATALOG')

                #Combine the HDUs into an HDUlist and save as a new file. Close the HDUlist.
                hdulist = fits.HDUList([prihdu, hdu_CATALOG, hdu_DELTA, hdu_VEL, hdu_COSMO])
                hdulist.writeto(filename,overwrite=True)
                hdulist.close

                success = 1

            elif output_format == 'picca':

                #Organise the data into picca-format arrays.
                picca_0 = self.DELTA_rows.T
                picca_1 = self.IVAR_rows.T
                picca_2 = self.LOGLAM_MAP

                picca_3_data = []
                for i in range(self.N_qso):
                    picca_3_data += [(self.RA[i],self.DEC[i],self.Z_QSO[i],self.PLATE[i],self.MJD[i],self.FIBER[i],self.THING_ID[i])]

                dtype = [('RA', '>f4'), ('DEC', '>f4'), ('Z', '>f4'), ('PLATE', '>f4'), ('MJD', '>f4'), ('FIBER', '>f4'), ('THING_ID', int)]
                picca_3 = np.array(picca_3_data,dtype=dtype)

                #Make the data into suitable HDUs.
                hdu_DELTA = fits.PrimaryHDU(data=picca_0,header=header)
                hdu_iv = fits.ImageHDU(data=picca_1,header=header,name='IV')
                hdu_LOGLAM_MAP = fits.ImageHDU(data=picca_2,header=header,name='LOGLAM_MAP')
                cols_CATALOG = fits.ColDefs(picca_3)
                hdu_CATALOG = fits.BinTableHDU.from_columns(cols_CATALOG,header=header,name='CATALOG')

                #Combine the HDUs into and HDUlist and save as a new file. Close the HDUlist.
                hdulist = fits.HDUList([hdu_DELTA, hdu_iv, hdu_LOGLAM_MAP, hdu_CATALOG])
                hdulist.writeto(filename,overwrite=True)
                hdulist.close()

                success = 1

            else:
                print('Output format "{}" not recognised.\nCurrent options are "colore" and "picca".'.format(output_format))
                output_format = raw_input('Please enter one of these options: ')

        return

    @classmethod
    def choose_qsos(cls,A,THING_IDs):

        rows = ['']*len(THING_IDs)
        s = set(THING_IDs)
        j = 0
        for i, qso in enumerate(A.THING_ID):
            if qso in s:
                rows[j] = i
                j=j+1

        N_qso = len(rows)
        N_cells = A.N_cells

        TYPE = A.TYPE[rows]
        RA = A.RA[rows]
        DEC = A.DEC[rows]
        Z_QSO = A.Z_QSO[rows]
        DZ_RSD = A.DZ_RSD[rows]
        THING_ID = A.THING_ID[rows]
        PLATE = A.PLATE[rows]
        MJD = A.MJD[rows]
        FIBER = A.FIBER[rows]

        DELTA_rows = A.DELTA_rows[rows,:]
        VEL_rows = A.VEL_rows[rows,:]
        IVAR_rows = A.IVAR_rows[rows,:]

        Z = A.Z
        LOGLAM_MAP = A.LOGLAM_MAP
        R = A.R
        D = A.D
        V = A.V

        return cls(N_qso,N_cells,TYPE,RA,DEC,Z_QSO,DZ_RSD,THING_ID,PLATE,MJD,FIBER,DELTA_rows,VEL_rows,IVAR_rows,R,Z,D,V,LOGLAM_MAP)

    #Maybe obsolete
    @classmethod
    def choose_cells(cls,A,cells):

        N_qso = A.N_qso
        N_cells = len(cells)

        TYPE = A.TYPE
        RA = A.RA
        DEC = A.DEC
        Z_QSO = A.Z_QSO
        DZ_RSD = A.DZ_RSD
        THING_ID = A.THING_ID
        PLATE = A.PLATE
        MJD = A.MJD
        FIBER = A.FIBER

        DELTA_rows = A.DELTA_rows[:,cells]
        VEL_rows = A.VEL_rows[:,cells]
        IVAR_rows = A.IVAR_rows[:,cells]

        Z = A.Z[cells]
        LOGLAM_MAP = A.LOGLAM_MAP[cells]
        R = A.R[cells]
        D = A.D[cells]
        V = A.V[cells]

        return cls(N_qso,N_cells,TYPE,RA,DEC,Z_QSO,DZ_RSD,THING_ID,PLATE,MJD,FIBER,DELTA_rows,VEL_rows,IVAR_rows,R,Z,D,V,LOGLAM_MAP)

    #Maybe obsolete
    @classmethod
    def crop(cls,A,THING_ID,cells):

        rows = ['']*len(THING_ID)
        s = set(THING_ID)
        j = 0
        for i, qso in enumerate(A.THING_ID):
            if qso in s:
                rows[j] = i
                j=j+1

        N_qso = len(rows)
        N_cells = len(cells)

        TYPE = A.TYPE[rows]
        RA = A.RA[rows]
        DEC = A.DEC[rows]
        Z_QSO = A.Z_QSO[rows]
        DZ_RSD = A.DZ_RSD[rows]
        THING_ID = A.THING_ID[rows]
        PLATE = A.PLATE[rows]
        MJD = A.MJD[rows]
        FIBER = A.FIBER[rows]

        DELTA_rows = A.DELTA_rows[rows,:]
        DELTA_rows = DELTA_rows[:,cells]

        VEL_rows = A.VEL_rows[rows,:]
        VEL_rows = VEL_rows[:,cells]

        IVAR_rows = A.IVAR_rows[rows,:]
        IVAR_rows = IVAR_rows[:,cells]

        Z = A.Z[cells]
        LOGLAM_MAP = A.LOGLAM_MAP[cells]
        R = A.R[cells]
        D = A.D[cells]
        V = A.V[cells]

        return cls(N_qso,N_cells,TYPE,RA,DEC,Z_QSO,DZ_RSD,THING_ID,PLATE,MJD,FIBER,DELTA_rows,VEL_rows,IVAR_rows,R,Z,D,V,LOGLAM_MAP)
