import numpy as np
from astropy.io import fits
from astropy.table import Table
import healpy as hp
import os
import argparse
import time

"""
To do:
 - add in options to trim pixels or make a random selection of quasars

 - inputs: input format, original_filename_structure, nodes, location, N_side, new_filename_structure, save_location, trim pixels?, reduce objects?, output format

 - it'd be nice to have labelled columns for node_pixel_map rather than relying on which number column we're looking at

 - need to add some useful outputs to use as progress checks

 - use keyword arguments for some things (output format, z_min, etc)

 - add options for merging pre-existing files:
    - never merge
    - merge if no duplicates, new file if duplicates (currently, this is implemented)
    - merge even if duplicates, obiously not creating duplicate entries
 - need something for the master file?

 - sub-function possibilities

 - change CoLoRe to colore to make things simpler?
"""


#Program to convert a .fits file of HEALPix pixels into several .fits files of smaller HEALPix pixels.
#Takes files of the format outputted by CoLoRe and converts into the same format that picca takes as input.

#Function to orchestrate the repixelisation process.
#First creates a master file recording the catalog of quasars being analysed.
#Then, creates per-HEALPix pixel files containing the density skewers in either colore or picca formats.

def split_pixels(original_location,original_filename_structure,file_numbers,input_format,N_side,z_min,save_location,new_filename_structure,output_format):

    #Determine number of nodes (1 per file).
    N_nodes = len(file_numbers)
    N_pixels = 12*N_side**2

    #Get list of nodes:
    nodes_int = file_numbers

    #Convert the numbers into strings of length 3. e.g. 14 -> '014'
    node_ID_length = 3
    nodes = numbers_to_strings(nodes_int,node_ID_length)

    #Set up master_ID_data.
    #This will become a list of tuples, with one tuple for each qso.
    #Each tuple contains the qso's pixel, and its THING_ID.
    master_ID_data = []

    #Set up the node_pixel_map.
    #node_pixel_map[i,j] is 1 if the jth node contains qsos from the ith pixel, 0 otherwise.
    node_pixel_map = np.zeros((N_pixels,max(nodes_int)+1))

    #Extract master_data from each file.
    for node in nodes:
        #Open the file and extract the angular coordinate data.
        h = fits.open(original_location + '/' + original_filename_structure.format(int(node)))
        N_qso = h[1].data.shape[0]

        if z_min > 0:
            indices = [i for i in range(N_qso) if h[1].data['Z_COSMO'][i]>z_min]
            RA = h[1].data['RA'][indices]
            DEC = h[1].data['DEC'][indices]
            N_qso = RA.shape[0]
        else:
            indices = list(range(N_qso))
            RA = h[1].data['RA']
            DEC = h[1].data['DEC']

        h.close()

        #Convert the coordinates to pixels and THING_IDs
        THING_ID = get_THING_ID(node,indices)
        pixel_ID = get_pixel_ID(N_side,RA,DEC)
        node_pixel_list = np.sort(list(set(pixel_ID)))

        #Add the THING_ID and pixel_ID from this node to the master list.
        for j in range(N_qso):
            master_ID_data += [(THING_ID[j],pixel_ID[j])]

        #Add information to node_pixel_map
        for pixel in node_pixel_list:
            if pixel >= 0:
                node_pixel_map[pixel,int(node)]=1


    #Sort the THING_IDs and pixel_IDs into the right order: first by pixel number, and then by THING_ID.
    dtype = [('THING_ID', 'U10'), ('PIX', int)]
    master_ID = np.array(master_ID_data, dtype=dtype)
    master_ID_sort = np.sort(master_ID, order=['PIX','THING_ID'])

    #Remove the quasars with invalid coordinates from the master file data.
    master_ID_sort_filter = master_ID_sort[master_ID_sort['PIX']>=0]

    #Make a list of the pixels that the files cover.
    pixel_list = list(set(master_ID_sort_filter['PIX']))
    N_pixels = len(pixel_list)

    print('Working on files named as:\n"{}",\nand numbered by nodes:\n{}\n'.format(original_filename_structure,nodes))
    print('\nOutput files will be named as:\n"{}",\nand numbered by zmin, nside and pixel:\n{}, {}, {}\n'.format(new_filename_structure,z_min,N_side,pixel_list))

    print('Creating a master file containing catalog data...')

    #Take the output from each splitting and make it into a master file
    master_header = fits.Header()
    master_header['NSIDE'] = N_side
    master = fits.BinTableHDU.from_columns(master_ID_sort_filter,header=master_header)

    #Make a primary HDU.
    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)

    #Make the .fits file.
    hdulist = fits.HDUList([prihdu,master])
    hdulist.writeto(save_location + '/' + 'zmin_{}_nside_{}_'.format(z_min,N_side) + 'master.fits')
    hdulist.close()

    print('Master file "{}" saved!\n'.format('zmin_{}_nside_{}_'.format(z_min,N_side) + 'master.fits'))

    #Make a file containing the THING_ID of all invalid coordinates.
    master_ID_sort_bad = master_ID_sort[master_ID_sort['PIX']<0]
    dtype = [('THING_ID', 'U10')]
    bad_THING_IDs = np.array(master_ID_sort_bad['THING_ID'],dtype=dtype)
    cols_bad_THING_IDs = fits.ColDefs(bad_THING_IDs)

    bad_master = fits.BinTableHDU.from_columns(cols_bad_THING_IDs,header=master_header)
    hdulist = fits.HDUList([prihdu,bad_master])
    hdulist.writeto(save_location + '/' + 'zmin_{}_nside_{}_'.format(z_min,N_side) + 'bad_coordinates.fits')
    hdulist.close()

    print('{} objects have been found to have invalid angular coordinates.'.format(bad_THING_IDs.shape[0]))
    print('THING_IDs have been stored in a file named {}.\n'.format('zmin_{}_nside_{}_'.format(z_min,N_side) + 'bad_coordinates.fits'))

    for pixel in pixel_list:

        print('Working on HEALPix pixel number {} ({} of {})...'.format(pixel,pixel_list.index(pixel)+1,len(pixel_list)))

        relevant_nodes = [node for node in nodes if node_pixel_map[pixel,int(node)]==1]

        nodes_included = 0
        combined_0 = []
        combined_1 = []
        combined_2 = []
        combined_3 = []
        combined_4 = []

        for node in relevant_nodes:

            print(' -> Extracting data from node {} ({} of {})...'.format(node,relevant_nodes.index(node)+1,len(relevant_nodes)))

            #Get the THING_IDs of the relevant (z_qso<z_min, in current pixel) quasars in the original file.
            relevant_THING_IDs = [qso['THING_ID'] for qso in master_ID_sort_filter if qso['PIX']==pixel and qso['THING_ID'][:3]==node]
            relevant_qsos = [qso for qso in master_ID_sort_filter if qso['PIX']==pixel and qso['THING_ID'][:3]==node]

            print('    -> {} relevant quasars found'.format(len(relevant_qsos)))

            #Convert this to the relevant rows of the quasars.
            relevant_rows = []
            for THING_ID in relevant_THING_IDs:
                relevant_rows += [int(THING_ID[-7:])]

            #Open the file corresponding to the current pixel and node.
            working_colore = fits.open(original_location + '/' + original_filename_structure.format(int(node)))

            #Trim any pixels corresponding to z<z_min.
            first_relevant_pixel = np.argmax(working_colore[4].data['Z'] >= z_min)

            #Extract the data from the file.
            working_colore_1 = working_colore[1].data[relevant_rows]
            working_colore_2 = working_colore[2].data[relevant_rows,first_relevant_pixel:]
            working_colore_3 = working_colore[3].data[relevant_rows,first_relevant_pixel:]
            working_colore_4 = working_colore[4].data[first_relevant_pixel:]

            #Close the file.
            working_colore.close()

            if output_format == 'picca':
                #Convert the colore data into the picca format.
                working_picca_0, working_picca_1, working_picca_2, working_picca_3 = colore_to_picca(relevant_THING_IDs,working_colore_1,working_colore_2,working_colore_3,working_colore_4)

                #Add the new data to the combined data.
                if nodes_included > 0:
                    combined_1, combined_2, combined_3, combined_4 = merge_picca(working_picca_0,working_picca_1,working_picca_2,working_picca_3,combined_0,combined_1,combined_2,combined_3)
                else:
                    combined_0 = working_picca_0
                    combined_1 = working_picca_1
                    combined_3 = working_picca_3
                    combined_2 = working_picca_2

            elif output_format == 'colore':
                #Add the new data to the combined data.
                if nodes_included > 0:
                    combined_1, combined_2, combined_3, combined_4 = merge_colore(working_colore_1,working_colore_2,working_colore_3,working_colore_4,combined_1,combined_2,combined_3,combined_4)
                else:
                    combined_1 = working_colore_1
                    combined_2 = working_colore_2
                    combined_3 = working_colore_3
                    combined_4 = working_colore_4

            nodes_included += 1

        if nodes_included > 0:

            #Set up some useful headers.
            header = fits.Header()
            header['NSIDE'] = N_side
            header['PIX'] = pixel
            header['LYA'] = 1215.67

            #See if file for this pixel already exists
            #If so, deal with it in a sensible way

            if output_format == 'picca':
                header['NQSO'] = combined_3.shape
                combined_filename = save_location + '/' + output_format + '_' + new_filename_structure.format(z_min,N_side,pixel)
                check_filename = combined_filename
                check = 0
                i = 2

                #Check to see if a file with this name already exists.
                while check==0:
                    if os.path.exists(check_filename):
                        print('File named {} already exists.'.format(check_filename))
                        #If so, open it and extract the data.
                        working_picca = fits.open(check_filename)
                        working_picca_0 = working_picca[0].data
                        working_picca_1 = working_picca[1].data
                        working_picca_2 = working_picca[2].data
                        working_picca_3 = working_picca[3].data
                        working_picca.close()

                        #Produce lists of the THING_IDs in the existing and new datasets.
                        working_THING_ID_list = list(working_picca_3['THING_ID'])
                        new_THING_ID_list = list(combined_3['THING_ID'])
                        combined_THING_ID_list = working_THING_ID_list + new_THING_ID_list
                        reduced_combined_THING_ID_list = list(set(list(combined_THING_ID_list)))

                        #Check for duplicates.
                        if len(reduced_combined_THING_ID_list) == len(working_THING_ID_list)+len(new_THING_ID_list):
                            #If none, then attempt to merge it with the current file.
                            print('No duplicate entries found. Data will be added to existing file.')
                            combined_1, combined_2, combined_3, combined_4 = merge_picca(working_picca_0,working_picca_1,working_picca_2,working_picca_3,combined_0,combined_1,combined_2,combined_3)
                            combined_filename = check_filename
                            check = 1
                        else:
                            check_filename = combined_filename[:-5] + '_' + str(i) + '.fits'
                            print('Duplicates entries found. Will attempt to save new data in a new file named {}.'.format(check_filename))
                            i = i+1
                    else:
                        combined_filename = check_filename
                        check = 1

                #Make the data into suitable HDUs.
                hdu_DELTA = fits.PrimaryHDU(data=combined_0,header=header)
                hdu_iv = fits.ImageHDU(data=combined_1,header=header,name='IV')
                hdu_LOGLAM_MAP = fits.ImageHDU(data=combined_2,header=header,name='LOGLAM_MAP')
                cols_CATALOG = fits.ColDefs(combined_3)
                hdu_CATALOG = fits.BinTableHDU.from_columns(cols_CATALOG,header=header,name='CATALOG')

                #Combine the HDUs into and HDUlist and save as a new file. Close the HDUlist.
                hdulist = fits.HDUList([hdu_DELTA, hdu_iv, hdu_LOGLAM_MAP, hdu_CATALOG])
                hdulist.writeto(combined_filename)

                hdulist.close()

                print('Finished file written!\n')


            elif output_format == 'colore':
                header['NQSO'] = combined_1.shape
                combined_filename = save_location + '/' + output_format + '_' + new_filename_structure.format(z_min,N_side,pixel)
                check_filename = combined_filename
                check = 0
                i = 2

                #Check to see if a file with this name already exists.
                while check==0:
                    if os.path.exists(check_filename):
                        print('File named {} already exists.'.format(check_filename))
                        #If so, open it and extract the data.
                        working_colore = fits.open(check_filename)
                        working_colore_1 = working_colore[1].data
                        working_colore_2 = working_colore[2].data
                        working_colore_3 = working_colore[3].data
                        working_colore_4 = working_colore[4].data
                        working_colore.close()

                        #Produce lists of the THING_IDs in the existing and new datasets.
                        working_THING_ID_list = list(working_colore_3['THING_ID'])
                        new_THING_ID_list = list(combined_3['THING_ID'])
                        combined_THING_ID_list = working_THING_ID_list + new_THING_ID_list
                        reduced_combined_THING_ID_list = list(set(list(combined_THING_ID_list)))

                        #Check for duplicates.
                        if len(reduced_combined_THING_ID_list) == len(working_THING_ID_list)+len(new_THING_ID_list):
                            #If none, then attempt to merge it with the current file.
                            print('No duplicate entries found. Data will be added to existing file.')
                            combined_1, combined_2, combined_3, combined_4 = merge_colore(working_colore_1,working_colore_2,working_colore_3,working_colore_4,combined_1,combined_2,combined_3,combined_4)
                            combined_filename = check_filename
                            check = 1
                        else:
                            print('Duplicates entries found. Will attempt to save new data in a new file named'.format(combined_filename))
                            check_filename = combined_filename[:-5] + '_' + str(i) + '.fits'
                            i = i+1
                    else:
                        combined_filename = check_filename
                        check = 1

                #Make the data into suitable HDUs.
                prihdr = fits.Header()
                prihdu = fits.PrimaryHDU(header=prihdr)
                cols_CATALOG = fits.ColDefs(combined_1)
                hdu_CATALOG = fits.BinTableHDU.from_columns(cols_CATALOG,header=header,name='CATALOG')
                hdu_DELTA = fits.ImageHDU(data=combined_2,header=header,name='DELTA')
                hdu_VEL = fits.ImageHDU(data=combined_3,header=header,name='VELOCITY')
                cols_COSMO = fits.ColDefs(combined_4)
                hdu_COSMO = fits.BinTableHDU(cols_COSMO,header=header,name='CATALOG')

                #Combine the HDUs into and HDUlist and save as a new file. Close the HDUlist.
                hdulist = fits.HDUList([prihdu, hdu_CATALOG, hdu_DELTA, hdu_VEL, hdu_COSMO])
                hdulist.writeto(combined_filename)

                hdulist.close()

                print('Finished file written!\n')

    return

#Function to produce a list of THING_IDs made of a 3-digit node number, then a 7-digit row number.
def get_THING_ID(node,row_numbers):

    N_qso = len(row_numbers)

    THING_ID = ['']*N_qso
    for i in range(N_qso):
        row_numbers[i] = str(row_numbers[i])
        if len(row_numbers[i])<=7:
            row_numbers[i] = '0'*(7-len(row_numbers[i]))+row_numbers[i]
        else:
            exit('The row number is too great to construct a unique THING_ID (more than 7 digits).')
        THING_ID[i] = node+row_numbers[i]

    return THING_ID

#Function to identify which HEALPix pixels a set of coordinates each live in for a given N_side.
def get_pixel_ID(N_side,RA,DEC):

    N_qso = RA.shape[0]

    #Convert DEC and RA in degrees to theta and phi in radians.
    theta = (np.pi/180.0)*(90.0-DEC)
    phi = (np.pi/180.0)*RA

    #Make a list of the HEALPix pixel coordinate of each quasar.
    pixel_ID = ['']*N_qso
    for i in range(N_qso):
        #Check that the angular coordinates are valid. Put all objects with invalid coordinates into a non-realistic ID number (-1).
        if 0 <= theta[i] <= np.pi and 0 <= phi[i] <= 2*np.pi:
            pixel_ID[i] = int(hp.pixelfunc.ang2pix(N_side,theta[i],phi[i]))
        else:
            pixel_ID[i] = -1

    return pixel_ID

#Function to merge two CoLoRe files.
def merge_colore(colore_1_1,colore_1_2,colore_1_3,colore_1_4,colore_2_1,colore_2_2,colore_2_3,colore_2_4):

    combined_1 = np.concatenate((colore_1_1,colore_2_1),axis=0)
    combined_2 = np.concatenate((colore_1_2,colore_2_2),axis=0)
    combined_3 = np.concatenate((colore_1_3,colore_2_3),axis=0)

    if colore_1_4 == colore_2_4:
        combined_4 = colore_1_4
    else:
        exit('Different background cosmologies, cannot combine files.')

    return combined_1, combined_2, combined_3, combined_4

#Function to merge two picca files.
def merge_picca(picca_1_0,picca_1_1,picca_1_2,picca_1_3,picca_2_0,picca_2_1,picca_2_2,picca_2_3):

    combined_0 = np.concatenate((picca_1_0,picca_2_0),axis=1)
    combined_1 = np.concatenate((picca_1_1,picca_2_1),axis=1)
    combined_3 = np.concatenate((picca_1_3,picca_2_3),axis=0)

    if picca_1_2 == picca_2_2:
        combined_2 = picca_1_2
    else:
        exit('Different loglam maps, cannot combine files.')

    return combined_0

#Function to convert a list of numbers to a list of n-digit strings.
def numbers_to_strings(numbers,string_length):

    strings=[]
    for number in numbers:
        number_str = str(number)
        if len(number_str)<=string_length:
            string = '0'*(string_length-len(number_str))+number_str
        else:
            exit('The node number is too great to construct a unique THING_ID (more than 3 digits).')

        strings += [string]

    return strings

#Function that takes the data from a CoLoRe format and a list of THING_IDs, and converts to the data in picca format.
def colore_to_picca(THING_IDs,colore_1,colore_2,colore_3,colore_4):

    #Read the information we need from the colore data arrays.
    RA = colore_1['RA']
    DEC = colore_1['DEC']
    z_qso = colore_1['Z_COSMO']
    N_qso = colore_2.shape[0]
    N_cells = colore_2.shape[1]
    z = colore_4['Z']

    #Add in extra information.
    #Can we get the lya line from a list of constants somewhere instead?
    lya = 1215.67
    iv = np.ones((N_qso,N_cells))
    PLATE = np.zeros(N_qso)
    MJD = np.zeros(N_qso)
    FIBER = np.zeros(N_qso)

    #Construct the picca data arrays.
    picca_0 = colore_2.T
    picca_1 = iv.T
    picca_2 = np.log10(lya*(1+z))

    picca_3_data = ['']*N_qso
    for i in range(N_qso):
        picca_3_data[i] = (RA[i],DEC[i],z_qso[i],PLATE[i],MJD[i],FIBER[i],THING_IDs[i])

    dtype = [('RA', '>f4'), ('DEC', '>f4'), ('Z', '>f4'), ('PLATE', '>f4'), ('MJD', '>f4'), ('FIBER', '>f4'), ('THING_ID', 'U10')]
    picca_3 = np.array(picca_3_data,dtype=dtype)

    return picca_0, picca_1, picca_2, picca_3
