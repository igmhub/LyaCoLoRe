import numpy as np
from astropy.io import fits
import healpy as hp
import os
import time

# TODO: remove SIGMA_G from the headers of the saved files as it cannot be relied upon.

lya = 1215.67

#Function to create a 'simulation_data' object given a specific pixel, information about the complete simulation, and the location/filenames of data files.
def make_gaussian_pixel_object(pixel,original_file_location,original_filename_structure,input_format,MOCKID_lookup,lambda_min=0,IVAR_cutoff=lya):

    #Determine which file numbers we need to look at for the current pixel.
    relevant_keys = [key for key in MOCKID_lookup.keys() if key[1]==pixel]
    files_included = 0

    #For each relevant file, extract the data and aggregate over all files into a 'combined' object.
    for key in relevant_keys:
        #Get the MOCKIDs of the relevant quasars: those that are located in the current pixel, stored in the current file.
        file_number = key[0]
        relevant_MOCKIDs = MOCKID_lookup[key]
        N_relevant_qso = len(relevant_MOCKIDs)

        #If there are some relevant quasars, open the data file and make it into a simulation_data object.
        #We use simulation_data.get_reduced_data to avoid loading all of the file's data into the object.
        if N_relevant_qso > 0:
            filename = original_file_location + '/' + original_filename_structure.format(file_number)
            working = simulation_data.get_gaussian_skewers(filename,file_number,input_format,relevant_MOCKIDs,lambda_min=lambda_min,IVAR_cutoff=IVAR_cutoff)

        #Combine the data from the working file with that from the files already looked at.
        if files_included > 0:
            combined = simulation_data.combine_files(combined,working,gaussian_only=True)
            files_included += 1
        else:
            combined = working
            files_included += 1

    pixel_object = combined

    return pixel_object

#Function to create a file structure based on a set of numbers, of the form "x//100"/"x".
def make_file_structure(base_location,numbers):

    first_level = []
    for number in numbers:
        first_level += [number//100]

    first_level_set = list(sorted(set(first_level)))

    for i in first_level_set:

        os.mkdir(base_location+str(i))

        for j, number in enumerate(numbers):

            if first_level[j] == i:
                os.mkdir(base_location+str(i)+'/'+str(number))

    return

#Function to convert a list of numbers to a list of n-digit strings.
def number_to_string(number,string_length):

    number_str = str(number)
    if len(number_str)<=string_length:
        string = '0'*(string_length-len(number_str))+number_str
    else:
        exit('The file number is too great to construct a unique MOCKID (more than 3 digits).')

    return string

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
            MOCKID = make_MOCKID(file_number,row_numbers)
    elif input_format == 'gaussian_colore':
        try:
            MOCKID = h[1].data['MOCKID']
        except KeyError:
            h_N_qso = h[1].data.shape[0]
            row_numbers = list(range(h_N_qso))
            MOCKID = make_MOCKID(file_number,row_numbers)
    elif input_format == 'picca':
        MOCKID = h[3].data['MOCKID']
    elif input_format == 'ID':
        MOCKID = h[1].data['MOCKID']

    return MOCKID

#Function to construct an array of MOCKIDs given a file_number and a list of row_numbers.
def make_MOCKID(file_number,row_numbers):

    N_qso = len(row_numbers)
    node = '0'*(len(str(file_number))-5) + str(file_number)

    MOCKID = ['']*N_qso
    for i in range(N_qso):
        row_numbers[i] = str(row_numbers[i])
        if len(row_numbers[i])<=7:
            row_numbers[i] = '0'*(7-len(row_numbers[i]))+row_numbers[i]
        else:
            exit('The row number is too great to construct a unique MOCKID (more than 7 digits).')
        MOCKID[i] = int(node+row_numbers[i])

    MOCKID = np.array(MOCKID)

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
def get_ID_data(original_file_location,original_filename_structure,file_number,input_format,N_side):

    ID_data = []
    cosmology = []
    N_pixels = 12*N_side**2

    #Open the file and extract the angular coordinate data.
    filename = original_file_location + '/' + original_filename_structure.format(file_number)
    h = fits.open(filename)

    #Extract the component parts of the master file's data from h.
    RA = get_RA(h,input_format)
    DEC = get_DEC(h,input_format)
    Z_QSO_NO_RSD = get_Z_QSO(h,input_format)
    DZ_RSD = get_DZ_RSD(h,input_format)
    MOCKID = get_MOCKID(h,input_format,file_number)
    h_R, h_Z, h_D, h_V = get_COSMO(h,input_format)

    h.close()

    #Construct the remaining component parts of the master file's data.
    pixel_ID = make_pixel_ID(N_side,RA,DEC)

    #Calculate Z_QSO_RSD.
    Z_QSO_RSD = Z_QSO_NO_RSD + DZ_RSD

    #Join the pieces of the ID_data together.
    ID_data = list(zip(RA,DEC,Z_QSO_NO_RSD,Z_QSO_RSD,MOCKID,pixel_ID))

    #Make file-pixel map element and MOCKID lookup.
    pixel_ID_set = list(sorted(set([pixel for pixel in pixel_ID if pixel>=0])))
    file_pixel_map_element = np.zeros(N_pixels)
    MOCKID_lookup_element = {}
    for pixel in pixel_ID_set:
        file_pixel_map_element[pixel] = 1
        MOCKID_pixel_list = [MOCKID[i] for i in range(len(pixel_ID)) if pixel_ID[i]==pixel]
        MOCKID_lookup_element = {**MOCKID_lookup_element,**{(file_number,pixel):MOCKID_pixel_list}}

    #Sort the MOCKIDs and pixel_IDs into the right order: first by pixel number, and then by MOCKID.
    dtype = [('RA', '>f4'), ('DEC', '>f4'), ('Z_QSO_NO_RSD', '>f4'), ('Z_QSO_RSD', '>f4'), ('MOCKID', int), ('PIXNUM', int)]
    ID = np.array(ID_data, dtype=dtype)
    ID_sort = np.sort(ID, order=['PIXNUM','MOCKID'])

    #Construct the cosmology array.
    cosmology_data = list(zip(h_R,h_Z,h_D,h_V))
    dtype = [('R', '>f4'), ('Z', '>f4'), ('D', '>f4'), ('V', '>f4')]
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
    dtype = [('RA','>f4'),('DEC','>f4'),('Z','>f4'),('THING_ID',int),('MJD','>f4'),('FIBERID',int),('PLATE',int),('PIXNUM',int)]
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

#From lya_mock_p1d.py
def get_tau(z,density,A,alpha=1.0):
    """transform lognormal density to optical depth, at each z"""
    # add redshift evolution to mean optical depth
    A = 0.374*pow((1+z)/4.0,5.10)

    TAU_rows = A*(density**alpha)

    return TAU_rows

#Function to make ivar mask
def make_IVAR_rows(lya,Z_QSO,LOGLAM_MAP):

    N_cells = LOGLAM_MAP.shape[0]
    N_qso = Z_QSO.shape[0]

    lya_lambdas = lya*(1+Z_QSO)
    IVAR_rows = np.ones((N_qso,N_cells),dtype='float32')
    lambdas = 10**LOGLAM_MAP

    for i in range(N_qso):

        last_relevant_cell = (np.argmax(lambdas > lya_lambdas[i]) - 1) % N_cells

        for j in range(last_relevant_cell+1,N_cells):
            IVAR_rows[i,j] = 0.

    return IVAR_rows

#Function to convert lognormal delta skewers (in rows) to gaussian field skewers (in rows).
def lognormal_to_gaussian(LN_DENSITY_DELTA_rows,SIGMA_G,D):

    LN_DENSITY_rows = 1.0 + LN_DENSITY_DELTA_rows

    GAUSSIAN_DELTA_rows = np.zeros(LN_DENSITY_DELTA_rows.shape)

    for j in range(GAUSSIAN_DELTA_rows.shape[1]):
        GAUSSIAN_DELTA_rows[:,j] = (np.log(LN_DENSITY_rows[:,j]))/D[j] + (D[j])*(SIGMA_G**2)/2

    GAUSSIAN_DELTA_rows = GAUSSIAN_DELTA_rows.astype('float32')

    return GAUSSIAN_DELTA_rows

#Function to convert gaussian field skewers (in rows) to lognormal delta skewers (in rows).
def gaussian_to_lognormal(GAUSSIAN_DELTA_rows,SIGMA_G,D):

    LN_DENSITY_rows = np.zeros(GAUSSIAN_DELTA_rows.shape)
    LN_DENSITY_DELTA_rows = np.zeros(GAUSSIAN_DELTA_rows.shape)

    for j in range(GAUSSIAN_DELTA_rows.shape[1]):
        LN_DENSITY_rows[:,j] = np.exp(D[j]*GAUSSIAN_DELTA_rows[:,j] - ((D[j])**2)*(SIGMA_G**2)/2.)

    LN_DENSITY_DELTA_rows = LN_DENSITY_rows - 1

    LN_DENSITY_DELTA_rows = LN_DENSITY_DELTA_rows.astype('float32')

    return LN_DENSITY_DELTA_rows

#Function to determine the first index corresponding to a value in an array greater than a minimum value.
def get_first_relevant_index(minimum,values):

    if minimum > 0:
        first_relevant_index = np.argmax(values >= minimum)
    else:
        first_relevant_index = 0

    return first_relevant_index

#Function to determine the indices corresponding to the values in an array greater than a minimum value.
def get_relevant_indices(minimum,values):

    N = values.shape[0]
    relevant_indices = [i for i in range(N) if values[i] > minimum]

    return relevant_indices

# TODO: sort out the float/integer problem
#Function to retrieve relevant simulation parameters from the param.cfg file.
def get_simulation_parameters(location,filename):

    #Create a string of the parameter file to search.
    parameter_str = open(location + '/' + filename,'r').read()

    #Define the parameters to search for and the intricacies of the parameter file.
    divider = '\n'
    parameters = [('dens_type','int'),('omega_M','>f4'),('omega_L','>f4'),('omega_B','>f4'),('h','>f4'),('w','>f4'),('ns','>f4'),('sigma_8','>f4')]
    equality_format = ' = '
    N_parameters = len(parameters)

    parameter_values = np.array([tuple([0]*N_parameters)],dtype=parameters)

    #For each parameter defined, find its value and put it into the output array.
    for i, parameter in enumerate(parameters):
        search_term = parameter[0] + equality_format
        data_type = parameter[1]
        N_char = len(search_term)

        first_char = parameter_str.find(search_term)
        parameter_values[0][i] = float(parameter_str[first_char+N_char:parameter_str.find(divider,first_char+N_char)])

    return parameter_values

#Function to normalise a set of delta skewer rows to zero mean according to given weights.
#If all weights for a given cell are zero, then the output will be zero in that cell for all skewers.
def normalise_deltas(DELTA_rows,weights):

    N_cells = DELTA_rows.shape[1]
    DELTA_rows_mean = np.zeros(N_cells)
    DELTA_rows_normalised = np.zeros(DELTA_rows.shape)

    for j in range(N_cells):
        if np.sum(weights[:,j]) != 0:
            DELTA_rows_mean[j] = np.average(DELTA_rows[:,j],weights=weights[:,j])
            DELTA_rows_normalised[:,j] = (DELTA_rows[:,j] + 1)/(DELTA_rows_mean[j] + 1) - 1

    return DELTA_rows_normalised

# TODO: write this
class simulation_parameters:
    def __init__():

        return

    @classmethod
    def get_parameters(cls,location,filename):

        return

#Function to calculate the mean of deltas, mean of deltas^2, and N.
def return_means(DELTA_rows,weights,sample_pc=1.0):
    DELTA_SQUARED_rows = DELTA_rows**2

    N = np.sum(weights)
    mean_DELTA = np.average(DELTA_rows,weights=weights)
    mean_DELTA_SQUARED = np.average(DELTA_SQUARED_rows,weights=weights)

    return N, mean_DELTA, mean_DELTA_SQUARED

#Function to take a list of sets of statistics (as produced by 'get_statistics'), and calculate means and variances.
def combine_means(means_list):

    means_shape = means_list[0].shape
    means_data_type = means_list[0].dtype

    quantities = ['GAUSSIAN_DELTA','GAUSSIAN_DELTA_SQUARED','DENSITY_DELTA','DENSITY_DELTA_SQUARED','F','F_SQUARED','F_DELTA','F_DELTA_SQUARED']

    combined_means = np.zeros(means_shape,dtype=means_data_type)
    for means_array in means_list:
        combined_means['N'] += means_array['N']

    for quantity in quantities:
        for means_array in means_list:
            combined_means[quantity] += (means_array['N']*means_array[quantity])/combined_means['N']

    return combined_means

#Function to convert a set of means of quantities and quantities squared (as outputted by 'combine_means') to a set of means and variances.
def means_to_statistics(means):

    statistics_dtype = [('N', '>f4')
        , ('GAUSSIAN_DELTA_MEAN', '>f4'), ('GAUSSIAN_DELTA_VAR', '>f4')
        , ('DENSITY_DELTA_MEAN', '>f4'), ('DENSITY_DELTA_VAR', '>f4')
        , ('F_MEAN', '>f4'), ('F_VAR', '>f4')
        , ('F_DELTA_MEAN', '>f4'), ('F_DELTA_VAR', '>f4')]

    statistics = np.zeros(means.shape,dtype=statistics_dtype)

    statistics['N'] = means['N']
    statistics['GAUSSIAN_DELTA_MEAN'] = means['GAUSSIAN_DELTA']
    statistics['DENSITY_DELTA_MEAN'] = means['DENSITY_DELTA']
    statistics['F_MEAN'] = means['F']
    statistics['F_DELTA_MEAN'] = means['F_DELTA']

    statistics['GAUSSIAN_DELTA_VAR'] = means['GAUSSIAN_DELTA_SQUARED'] - means['GAUSSIAN_DELTA']**2
    statistics['DENSITY_DELTA_VAR'] = means['DENSITY_DELTA_SQUARED'] - means['DENSITY_DELTA']**2
    statistics['F_VAR'] = means['F_SQUARED'] - means['F']**2
    statistics['F_DELTA_VAR'] = means['F_DELTA_SQUARED'] - means['F_DELTA']**2

    return statistics

#Function to write the statistics data to file, along with an HDU extension contanint cosmology data.
def write_statistics(location,N_side,statistics,cosmology_data):

    #Construct HDU from the statistics array.
    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)
    cols_stats = fits.ColDefs(statistics)
    hdu_stats = fits.BinTableHDU.from_columns(cols_stats,name='STATISTICS')
    cols_cosmology = fits.ColDefs(cosmology_data)
    hdu_cosmology = fits.BinTableHDU.from_columns(cols_cosmology,name='COSMO')

    #Put the HDU into an HDUlist and save as a new file. Close the HDUlist.
    filename = 'nside_{}_statistics.fits'.format(N_side)
    hdulist = fits.HDUList([prihdu,hdu_stats,hdu_cosmology])
    hdulist.writeto(location+filename)
    hdulist.close

    return

#Definition of a generic 'simulation_data' class, from which it is easy to save in new formats.
class simulation_data:
    #Initialisation function.
    def __init__(self,N_qso,N_cells,SIGMA_G,ALPHA,TYPE,RA,DEC,Z_QSO,DZ_RSD,MOCKID,PLATE,MJD,FIBER,GAUSSIAN_DELTA_rows,DENSITY_DELTA_rows,VEL_rows,IVAR_rows,F_rows,R,Z,D,V,LOGLAM_MAP,A):

        self.N_qso = N_qso
        self.N_cells = N_cells
        self.SIGMA_G = SIGMA_G
        self.ALPHA = ALPHA

        self.TYPE = TYPE
        self.RA = RA
        self.DEC = DEC
        self.Z_QSO = Z_QSO
        self.DZ_RSD = DZ_RSD
        self.MOCKID = MOCKID
        self.PLATE = PLATE
        self.MJD = MJD
        self.FIBER = FIBER

        self.GAUSSIAN_DELTA_rows = GAUSSIAN_DELTA_rows
        self.DENSITY_DELTA_rows = DENSITY_DELTA_rows
        self.VEL_rows = VEL_rows
        self.IVAR_rows = IVAR_rows
        self.F_rows = F_rows

        self.R = R
        self.Z = Z
        self.D = D
        self.V = V
        self.LOGLAM_MAP = LOGLAM_MAP
        self.A = A

        return

    #Method to extract reduced data from an input file of a given format, with a given list of MOCKIDs.
    @classmethod
    def get_gaussian_skewers(cls,filename,file_number,input_format,MOCKIDs=None,lambda_min=0,IVAR_cutoff=lya,SIGMA_G=None):

        lya = 1215.67

        h = fits.open(filename)

        h_MOCKID = get_MOCKID(h,input_format,file_number)
        h_R, h_Z, h_D, h_V = get_COSMO(h,input_format)
        h_lya_lambdas = get_lya_lambdas(h,input_format)

        if MOCKIDs != None:
            #Work out which rows in the hdulist we are interested in.
            rows = ['']*len(MOCKIDs)
            s = set(MOCKIDs)
            j = 0
            for i, qso in enumerate(h_MOCKID):
                if qso in s:
                    rows[j] = i
                    j = j+1
        else:
            rows = list(range(h_MOCKID.shape[0]))

        #Calculate the first_relevant_cell.
        first_relevant_cell = np.argmax(h_lya_lambdas >= lambda_min)
        actual_lambda_min = h_lya_lambdas[first_relevant_cell]

        if input_format == 'physical_colore':

            #Extract data from the HDUlist.
            TYPE = h[1].data['TYPE'][rows]
            RA = h[1].data['RA'][rows]
            DEC = h[1].data['DEC'][rows]
            Z_QSO = h[1].data['Z_COSMO'][rows]
            DZ_RSD = h[1].data['DZ_RSD'][rows]

            DENSITY_DELTA_rows = h[2].data[rows,first_relevant_cell:]

            VEL_rows = h[3].data[rows,first_relevant_cell:]

            Z = h[4].data['Z'][first_relevant_cell:]
            R = h[4].data['R'][first_relevant_cell:]
            D = h[4].data['D'][first_relevant_cell:]
            V = h[4].data['V'][first_relevant_cell:]

            #Derive the number of quasars and cells in the file.
            N_qso = RA.shape[0]
            N_cells = Z.shape[0]
            if SIGMA_G == None:
                SIGMA_G = h[4].header['SIGMA_G']

            #Derive the MOCKID and LOGLAM_MAP.
            MOCKID = get_MOCKID(h,input_format,file_number)
            LOGLAM_MAP = np.log10(lya*(1+Z))

            #Calculate the Gaussian skewers.
            GAUSSIAN_DELTA_rows = lognormal_to_gaussian(DENSITY_DELTA_rows,SIGMA_G,D)

            #Set the remaining variables to None
            A = None
            ALPHA = None
            TAU_rows = None
            F_rows = None
            DENSITY_DELTA_rows = None

            #Insert placeholder values for remaining variables.
            PLATE = MOCKID
            MJD = np.zeros(N_qso)
            FIBER = np.zeros(N_qso)

            IVAR_rows = make_IVAR_rows(IVAR_cutoff,Z_QSO,LOGLAM_MAP)

        elif input_format == 'gaussian_colore':

            #Extract data from the HDUlist.
            TYPE = h[1].data['TYPE'][rows]
            RA = h[1].data['RA'][rows]
            DEC = h[1].data['DEC'][rows]
            Z_QSO = h[1].data['Z_COSMO'][rows]
            DZ_RSD = h[1].data['DZ_RSD'][rows]

            GAUSSIAN_DELTA_rows = h[2].data[rows,first_relevant_cell:]

            VEL_rows = h[3].data[rows,first_relevant_cell:]

            Z = h[4].data['Z'][first_relevant_cell:]
            R = h[4].data['R'][first_relevant_cell:]
            D = h[4].data['D'][first_relevant_cell:]
            V = h[4].data['V'][first_relevant_cell:]

            #Derive the number of quasars and cells in the file.
            N_qso = RA.shape[0]
            N_cells = Z.shape[0]
            if SIGMA_G == None:
                SIGMA_G = h[4].header['SIGMA_G']

            #Derive the MOCKID and LOGLAM_MAP.
            MOCKID = get_MOCKID(h,input_format,file_number)
            LOGLAM_MAP = np.log10(lya*(1+Z))

            #Set the remaining variables to None
            A = None
            ALPHA = None
            TAU_rows = None
            F_rows = None
            DENSITY_DELTA_rows = None

            #Insert placeholder values for remaining variables.
            PLATE = MOCKID
            MJD = np.zeros(N_qso)
            FIBER = np.zeros(N_qso)

            IVAR_rows = make_IVAR_rows(IVAR_cutoff,Z_QSO,LOGLAM_MAP)

        elif input_format == 'picca':

            #Extract data from the HDUlist.
            DENSITY_DELTA_rows = h[0].data.T[rows,first_relevant_cell:]

            IVAR_rows = h[1].data.T[rows,first_relevant_cell:]

            LOGLAM_MAP = h[2].data[first_relevant_cell:]

            RA = h[3].data['RA'][rows]
            DEC = h[3].data['DEC'][rows]
            Z_QSO = h[3].data['Z'][rows]
            PLATE = h[3].data['PLATE'][rows]
            MJD = h[3].data['MJD'][rows]
            FIBER = h[3].data['FIBER'][rows]
            MOCKID = h[3].data['THING_ID'][rows]

            #Derive the number of quasars and cells in the file.
            N_qso = RA.shape[0]
            N_cells = LOGLAM_MAP.shape[0]
            if SIGMA_G == None:
                SIGMA_G = h[4].header['SIGMA_G']

            #Derive Z and transmitted flux fraction.
            Z = (10**LOGLAM_MAP)/lya - 1

            #Calculate the Gaussian skewers.
            GAUSSIAN_DELTA_rows = lognormal_to_gaussian(DENSITY_DELTA_rows,SIGMA_G,D)

            #Set the remaining variables to None
            A = None
            ALPHA = None
            TAU_rows = None
            F_rows = None
            DENSITY_DELTA_rows = None

            """
            Can we calculate DZ_RSD,R,D,V?
            """

            #Insert placeholder variables for remaining variables.
            TYPE = np.zeros(RA.shape[0])
            R = np.zeros(Z.shape[0])
            D = np.zeros(Z.shape[0])
            V = np.zeros(Z.shape[0])
            DZ_RSD = np.zeros(RA.shape[0])
            VEL_rows = np.zeros(DENSITY_DELTA_rows.shape)

        else:
            print('Input format not recognised: current options are "colore" and "picca".')
            print('Please choose one of these options and try again.')

        h.close()

        return cls(N_qso,N_cells,SIGMA_G,ALPHA,TYPE,RA,DEC,Z_QSO,DZ_RSD,MOCKID,PLATE,MJD,FIBER,GAUSSIAN_DELTA_rows,DENSITY_DELTA_rows,VEL_rows,IVAR_rows,F_rows,R,Z,D,V,LOGLAM_MAP,A)

    def add_small_scale_gaussian_fluctuations(self,dl,amplitude,white_noise=False):

        #Add small scale fluctuations
        #Redefine the necessary variables (N_cells, Z, D etc)
        #Warning if there are already physical/flux skewers

        return

    #Function to add physical skewers to the object via a lognormal transformation.
    def compute_physical_skewers(self,density_type='lognormal'):

        self.DENSITY_DELTA_rows = gaussian_to_lognormal(self.GAUSSIAN_DELTA_rows,self.SIGMA_G,self.D)

        return

    #Function to add flux skewers to the object.
    def compute_flux_skewers(self,A=None,alpha=1.0):

        if not A:
            A = 0.374*pow((1+self.Z)/4.0,5.10)

        self.TAU_rows = get_tau(self.Z,self.DENSITY_DELTA_rows+1,A,alpha=1.0)
        self.F_rows = np.exp(-self.TAU_rows)

        return

    #Method to combine data from two objects into one.
    # TODO: add something to check that we can just take values from 1 of the objects
    @classmethod
    def combine_files(cls,object_A,object_B,gaussian_only=False):

        N_qso = object_A.N_qso + object_B.N_qso

        """
        something to check N_cells is the same in both files
        """

        N_cells = object_A.N_cells
        SIGMA_G = object_A.SIGMA_G
        ALPHA = object_A.ALPHA

        TYPE = np.concatenate((object_A.TYPE,object_B.TYPE),axis=0)
        RA = np.concatenate((object_A.RA,object_B.RA),axis=0)
        DEC = np.concatenate((object_A.DEC,object_B.DEC),axis=0)
        Z_QSO = np.concatenate((object_A.Z_QSO,object_B.Z_QSO),axis=0)
        DZ_RSD = np.concatenate((object_A.DZ_RSD,object_B.DZ_RSD),axis=0)
        MOCKID = np.concatenate((object_A.MOCKID,object_B.MOCKID),axis=0)
        PLATE = np.concatenate((object_A.PLATE,object_B.PLATE),axis=0)
        MJD = np.concatenate((object_A.MJD,object_B.MJD),axis=0)
        FIBER = np.concatenate((object_A.FIBER,object_B.FIBER),axis=0)

        if not gaussian_only:
            GAUSSIAN_DELTA_rows = np.concatenate((object_A.GAUSSIAN_DELTA_rows,object_B.GAUSSIAN_DELTA_rows),axis=0)
            DENSITY_DELTA_rows = np.concatenate((object_A.DENSITY_DELTA_rows,object_B.DENSITY_DELTA_rows),axis=0)
            VEL_rows = np.concatenate((object_A.VEL_rows,object_B.VEL_rows),axis=0)
            IVAR_rows = np.concatenate((object_A.IVAR_rows,object_B.IVAR_rows),axis=0)
            F_rows = np.concatenate((object_A.F_rows,object_B.F_rows),axis=0)
        else:
            GAUSSIAN_DELTA_rows = np.concatenate((object_A.GAUSSIAN_DELTA_rows,object_B.GAUSSIAN_DELTA_rows),axis=0)
            DENSITY_DELTA_rows = None
            VEL_rows = np.concatenate((object_A.VEL_rows,object_B.VEL_rows),axis=0)
            IVAR_rows = np.concatenate((object_A.IVAR_rows,object_B.IVAR_rows),axis=0)
            F_rows = None

        """
        Something to check this is ok?
        """

        Z = object_A.Z
        LOGLAM_MAP = object_A.LOGLAM_MAP
        R = object_A.R
        D = object_A.D
        V = object_A.V
        A = object_A.A

        return cls(N_qso,N_cells,SIGMA_G,ALPHA,TYPE,RA,DEC,Z_QSO,DZ_RSD,MOCKID,PLATE,MJD,FIBER,GAUSSIAN_DELTA_rows,DENSITY_DELTA_rows,VEL_rows,IVAR_rows,F_rows,R,Z,D,V,LOGLAM_MAP,A)

    #Function to save data as a Gaussian colore file.
    def save_as_gaussian_colore(self,location,filename,header):

        #Organise the data into colore-format arrays.
        colore_1_data = []
        for i in range(self.N_qso):
            colore_1_data += [(self.TYPE[i],self.RA[i],self.DEC[i],self.Z_QSO[i],self.DZ_RSD[i],self.MOCKID[i])]

        dtype = [('TYPE', '>f4'), ('RA', '>f4'), ('DEC', '>f4'), ('Z_COSMO', '>f4'), ('DZ_RSD', '>f4'), ('MOCKID', int)]
        colore_1 = np.array(colore_1_data,dtype=dtype)
        colore_2 = self.GAUSSIAN_DELTA_rows
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
        hdu_GAUSSIAN = fits.ImageHDU(data=colore_2,header=header,name='GAUSSIAN_DELTA')
        hdu_VEL = fits.ImageHDU(data=colore_3,header=header,name='VELOCITY')
        cols_COSMO = fits.ColDefs(colore_4)
        hdu_COSMO = fits.BinTableHDU.from_columns(cols_COSMO,header=header,name='COSMO')

        #Combine the HDUs into an HDUlist and save as a new file. Close the HDUlist.
        hdulist = fits.HDUList([prihdu, hdu_CATALOG, hdu_GAUSSIAN, hdu_VEL, hdu_COSMO])
        hdulist.writeto(location+filename)
        hdulist.close

        return

    #Function to save data as a picca density file.
    def save_as_picca_gaussian(self,location,filename,header,zero_mean_delta=False,lambda_min=0):

        z_min = max(lambda_min/lya - 1, 0)
        lya_lambdas = 10**self.LOGLAM_MAP

        #Determine the first cell which corresponds to a lya_line at wavelength > lambda_min
        first_relevant_cell = get_first_relevant_index(lambda_min,lya_lambdas)

        #Determine the furthest cell which is still relevant: i.e. the one in which at least one QSO has non-zero value of IVAR.
        furthest_QSO_index = np.argmax(self.Z_QSO)
        #last_relevant_cell = (np.argmax(self.IVAR_rows[furthest_QSO_index,:]==0) - 1) % self.N_cells
        last_relevant_cell = self.N_cells - 1

        #Determine the relevant QSOs: those that have relevant cells (IVAR > 0) beyond the first_relevant_cell.
        #We impose a minimum number of cells per skewer here to avoid problems with picca.
        min_number_cells = 2
        relevant_QSOs = []
        for i in range(self.N_qso):
            min_relevant_cells = self.IVAR_rows[i,first_relevant_cell:(first_relevant_cell+min_number_cells)]
            if np.sum(min_relevant_cells)==min_number_cells:
                relevant_QSOs += [i]

        #Trim data according to the relevant cells and QSOs.
        relevant_GAUSSIAN_DELTA_rows = self.GAUSSIAN_DELTA_rows[relevant_QSOs,first_relevant_cell:last_relevant_cell+1]
        relevant_IVAR_rows = self.IVAR_rows[relevant_QSOs,first_relevant_cell:last_relevant_cell+1]
        relevant_LOGLAM_MAP = self.LOGLAM_MAP[first_relevant_cell:last_relevant_cell+1]

        #If desired, enforce that the Delta rows have zero mean.
        if zero_mean_delta == True:
            relevant_GAUSSIAN_DELTA_rows = normalise_deltas(relevant_GAUSSIAN_DELTA_rows,relevant_IVAR_rows)

        #Organise the data into picca-format arrays.
        picca_0 = relevant_GAUSSIAN_DELTA_rows.T
        picca_1 = relevant_IVAR_rows.T
        picca_2 = relevant_LOGLAM_MAP

        picca_3_data = []
        for i in range(self.N_qso):
            if i in relevant_QSOs:
                picca_3_data += [(self.RA[i],self.DEC[i],self.Z_QSO[i],self.PLATE[i],self.MJD[i],self.FIBER[i],self.MOCKID[i])]

        dtype = [('RA', '>f4'), ('DEC', '>f4'), ('Z', '>f4'), ('PLATE', int), ('MJD', '>f4'), ('FIBER', int), ('THING_ID', int)]
        picca_3 = np.array(picca_3_data,dtype=dtype)

        #Make the data into suitable HDUs.
        hdu_DELTA = fits.PrimaryHDU(data=picca_0,header=header)
        hdu_iv = fits.ImageHDU(data=picca_1,header=header,name='IV')
        hdu_LOGLAM_MAP = fits.ImageHDU(data=picca_2,header=header,name='LOGLAM_MAP')
        cols_CATALOG = fits.ColDefs(picca_3)
        hdu_CATALOG = fits.BinTableHDU.from_columns(cols_CATALOG,header=header,name='CATALOG')

        #Combine the HDUs into and HDUlist and save as a new file. Close the HDUlist.
        hdulist = fits.HDUList([hdu_DELTA, hdu_iv, hdu_LOGLAM_MAP, hdu_CATALOG])
        hdulist.writeto(location+filename)
        hdulist.close()

        return

    #Function to save data as a Lognormal colore file.
    def save_as_physical_colore(self,location,filename,header):

        #Organise the data into colore-format arrays.
        colore_1_data = []
        for i in range(self.N_qso):
            colore_1_data += [(self.TYPE[i],self.RA[i],self.DEC[i],self.Z_QSO[i],self.DZ_RSD[i],self.MOCKID[i])]

        dtype = [('TYPE', '>f4'), ('RA', '>f4'), ('DEC', '>f4'), ('Z_COSMO', '>f4'), ('DZ_RSD', '>f4'), ('MOCKID', int)]
        colore_1 = np.array(colore_1_data,dtype=dtype)

        colore_2 = self.DENSITY_DELTA_rows
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
        hdu_DELTA = fits.ImageHDU(data=colore_2,header=header,name='PHYSICAL_DELTA')
        hdu_VEL = fits.ImageHDU(data=colore_3,header=header,name='VELOCITY')
        cols_COSMO = fits.ColDefs(colore_4)
        hdu_COSMO = fits.BinTableHDU.from_columns(cols_COSMO,header=header,name='COSMO')

        #Combine the HDUs into an HDUlist and save as a new file. Close the HDUlist.
        hdulist = fits.HDUList([prihdu, hdu_CATALOG, hdu_DELTA, hdu_VEL, hdu_COSMO])
        hdulist.writeto(location+filename)
        hdulist.close

        return

    #Function to save data as a picca density file.
    def save_as_picca_density(self,location,filename,header,zero_mean_delta=False,lambda_min=0):

        z_min = max(lambda_min/lya - 1, 0)
        lya_lambdas = 10**self.LOGLAM_MAP

        #Determine the first cell which corresponds to a lya_line at wavelength > lambda_min
        first_relevant_cell = get_first_relevant_index(lambda_min,lya_lambdas)

        #Determine the furthest cell which is still relevant: i.e. the one in which at least one QSO has non-zero value of IVAR.
        furthest_QSO_index = np.argmax(self.Z_QSO)
        #last_relevant_cell = (np.argmax(self.IVAR_rows[furthest_QSO_index,:]==0) - 1) % self.N_cells
        last_relevant_cell = self.N_cells - 1

        #Determine the relevant QSOs: those that have relevant cells (IVAR > 0) beyond the first_relevant_cell.
        #We impose a minimum number of cells per skewer here to avoid problems with picca.
        min_number_cells = 2
        relevant_QSOs = []
        for i in range(self.N_qso):
            min_relevant_cells = self.IVAR_rows[i,first_relevant_cell:(first_relevant_cell+min_number_cells)]
            if np.sum(min_relevant_cells)==min_number_cells:
                relevant_QSOs += [i]

        #Trim data according to the relevant cells and QSOs.
        relevant_DENSITY_DELTA_rows = self.DENSITY_DELTA_rows[relevant_QSOs,first_relevant_cell:last_relevant_cell+1]
        relevant_IVAR_rows = self.IVAR_rows[relevant_QSOs,first_relevant_cell:last_relevant_cell+1]
        relevant_LOGLAM_MAP = self.LOGLAM_MAP[first_relevant_cell:last_relevant_cell+1]

        #If desired, enforce that the Delta rows have zero mean.
        if zero_mean_delta == True:
            relevant_DENSITY_DELTA_rows = normalise_deltas(relevant_DENSITY_DELTA_rows,relevant_IVAR_rows)

        #Organise the data into picca-format arrays.
        picca_0 = relevant_DENSITY_DELTA_rows.T
        picca_1 = relevant_IVAR_rows.T
        picca_2 = relevant_LOGLAM_MAP

        picca_3_data = []
        for i in range(self.N_qso):
            if i in relevant_QSOs:
                picca_3_data += [(self.RA[i],self.DEC[i],self.Z_QSO[i],self.PLATE[i],self.MJD[i],self.FIBER[i],self.MOCKID[i])]

        dtype = [('RA', '>f4'), ('DEC', '>f4'), ('Z', '>f4'), ('PLATE', int), ('MJD', '>f4'), ('FIBER', int), ('THING_ID', int)]
        picca_3 = np.array(picca_3_data,dtype=dtype)

        #Make the data into suitable HDUs.
        hdu_DELTA = fits.PrimaryHDU(data=picca_0,header=header)
        hdu_iv = fits.ImageHDU(data=picca_1,header=header,name='IV')
        hdu_LOGLAM_MAP = fits.ImageHDU(data=picca_2,header=header,name='LOGLAM_MAP')
        cols_CATALOG = fits.ColDefs(picca_3)
        hdu_CATALOG = fits.BinTableHDU.from_columns(cols_CATALOG,header=header,name='CATALOG')

        #Combine the HDUs into and HDUlist and save as a new file. Close the HDUlist.
        hdulist = fits.HDUList([hdu_DELTA, hdu_iv, hdu_LOGLAM_MAP, hdu_CATALOG])
        hdulist.writeto(location+filename)
        hdulist.close()

        return

    #Function to save data as a transmission file.
    def save_as_transmission(self,location,filename,header,lambda_min=0):

        z_min = max(lambda_min/lya - 1, 0)
        lya_lambdas = 10**self.LOGLAM_MAP

        first_relevant_cell = get_first_relevant_index(lambda_min,lya_lambdas)

        #Update lambda_min and z_min to the values of lambda and z of the first relevant cell.
        #This avoids including quasars that have Z_QSO higher than the original z_min, but no relevant skewer cells.
        lambda_min = lya_lambdas[first_relevant_cell]
        z_min = self.Z[first_relevant_cell]

        relevant_QSOs = get_relevant_indices(z_min,self.Z_QSO)

        #transmission_1_data = []
        #for i in relevant_QSOs:
        #    transmission_1_data += [(self.RA[i],self.DEC[i],self.Z_QSO[i],self.MOCKID[i])]

        transmission_1_data = list(zip(self.RA[relevant_QSOs],self.DEC[relevant_QSOs],self.Z_QSO[relevant_QSOs],self.MOCKID[relevant_QSOs]))

        dtype = [('RA', '>f4'), ('DEC', '>f4'), ('Z', '>f4'), ('MOCKID', int)]
        transmission_1 = np.array(transmission_1_data,dtype=dtype)

        transmission_2 = 10**(self.LOGLAM_MAP)[first_relevant_cell:]
        transmission_3 = self.F_rows[relevant_QSOs,first_relevant_cell:].T

        #Construct HDUs from the data arrays.
        prihdr = fits.Header()
        prihdu = fits.PrimaryHDU(header=prihdr)
        cols_METADATA = fits.ColDefs(transmission_1)
        hdu_METADATA = fits.BinTableHDU.from_columns(cols_METADATA,header=header,name='METADATA')
        hdu_WAVELENGTH = fits.ImageHDU(data=transmission_2,header=header,name='WAVELENGTH')
        hdu_TRANSMISSION = fits.ImageHDU(data=transmission_3,header=header,name='TRANSMISSION')

        #Combine the HDUs into and HDUlist and save as a new file. Close the HDUlist.
        hdulist = fits.HDUList([prihdu, hdu_METADATA, hdu_WAVELENGTH, hdu_TRANSMISSION])
        hdulist.writeto(location+filename)
        hdulist.close()

        return

    #Function to save data as a picca flux file.
    def save_as_picca_flux(self,location,filename,header,lambda_min=0):

        z_min = max(lambda_min/lya - 1, 0)
        lya_lambdas = 10**self.LOGLAM_MAP

        #Determine the first cell which corresponds to a lya_line at wavelength > lambda_min
        first_relevant_cell = get_first_relevant_index(lambda_min,lya_lambdas)

        #Determine the furthest cell which is still relevant: i.e. the one in which at least one QSO has non-zero value of IVAR.
        furthest_QSO_index = np.argmax(self.Z_QSO)
        #last_relevant_cell = (np.argmax(self.IVAR_rows[furthest_QSO_index,:]==0) - 1) % self.N_cells
        last_relevant_cell = self.N_cells - 1

        #Determine the relevant QSOs: those that have relevant cells (IVAR > 0) beyond the first_relevant_cell.
        #We impose a minimum number of cells per skewer here to avoid problems with picca.
        min_number_cells = 2
        relevant_QSOs = []
        for i in range(self.N_qso):
            min_relevant_cells = self.IVAR_rows[i,first_relevant_cell:(first_relevant_cell+min_number_cells)]
            if np.sum(min_relevant_cells)==min_number_cells:
                relevant_QSOs += [i]

        #Trim data according to the relevant cells and QSOs.
        relevant_F_rows = self.F_rows[relevant_QSOs,first_relevant_cell:last_relevant_cell+1]
        relevant_IVAR_rows = self.IVAR_rows[relevant_QSOs,first_relevant_cell:last_relevant_cell+1]
        relevant_LOGLAM_MAP = self.LOGLAM_MAP[first_relevant_cell:last_relevant_cell+1]

        #Calculate mean F as a function of z for the relevant cells, then F_DELTA_rows.
        #This is done with a 'hack' to avoid problems with weights summing to zero.
        # TODO: find a neater way to deal with this
        small = 1.0e-10
        relevant_F_BAR = np.average(relevant_F_rows,weights=relevant_IVAR_rows+small,axis=0)
        relevant_F_DELTA_rows = ((relevant_F_rows)/relevant_F_BAR - 1)*relevant_IVAR_rows

        #Organise the data into picca-format arrays.
        picca_0 = relevant_F_DELTA_rows.T
        picca_1 = relevant_IVAR_rows.T
        picca_2 = relevant_LOGLAM_MAP

        picca_3_data = []
        for i in range(self.N_qso):
            if i in relevant_QSOs:
                picca_3_data += [(self.RA[i],self.DEC[i],self.Z_QSO[i],self.PLATE[i],self.MJD[i],self.FIBER[i],self.MOCKID[i])]

        dtype = [('RA', '>f4'), ('DEC', '>f4'), ('Z', '>f4'), ('PLATE', int), ('MJD', '>f4'), ('FIBER', int), ('THING_ID', int)]
        picca_3 = np.array(picca_3_data,dtype=dtype)

        #Make the data into suitable HDUs.
        hdu_F = fits.PrimaryHDU(data=picca_0,header=header)
        hdu_iv = fits.ImageHDU(data=picca_1,header=header,name='IV')
        hdu_LOGLAM_MAP = fits.ImageHDU(data=picca_2,header=header,name='LOGLAM_MAP')
        cols_CATALOG = fits.ColDefs(picca_3)
        hdu_CATALOG = fits.BinTableHDU.from_columns(cols_CATALOG,header=header,name='CATALOG')

        #Combine the HDUs into and HDUlist and save as a new file. Close the HDUlist.
        hdulist = fits.HDUList([hdu_F, hdu_iv, hdu_LOGLAM_MAP, hdu_CATALOG])
        hdulist.writeto(location+filename)
        hdulist.close()

        return

    #Function to save the mean and variance of the different quantities as a function of Z.
    def get_means(self,lambda_min=0):

        #Determine the relevant cells and QSOs.
        lya_lambdas = 10**self.LOGLAM_MAP

        #Determine the first cell which corresponds to a lya_line at wavelength > lambda_min
        first_relevant_cell = get_first_relevant_index(lambda_min,lya_lambdas)

        #Determine the furthest cell which is still relevant: i.e. the one in which at least one QSO has non-zero value of IVAR.
        furthest_QSO_index = np.argmax(self.Z_QSO)
        #last_relevant_cell = (np.argmax(self.IVAR_rows[furthest_QSO_index,:]==0) - 1) % self.N_cells
        last_relevant_cell = self.N_cells - 1

        #Determine the relevant QSOs: those that have relevant cells (IVAR > 0) beyond the first_relevant_cell.
        relevant_QSOs = [i for i in range(self.N_qso) if self.IVAR_rows[i,first_relevant_cell] == 1]

        #Trim data according to the relevant cells and QSOs.
        relevant_DENSITY_DELTA_rows = self.DENSITY_DELTA_rows[relevant_QSOs,first_relevant_cell:last_relevant_cell+1]
        relevant_GAUSSIAN_DELTA_rows = self.GAUSSIAN_DELTA_rows[relevant_QSOs,first_relevant_cell:last_relevant_cell+1]
        relevant_F_rows = self.F_rows[relevant_QSOs,first_relevant_cell:last_relevant_cell+1]
        relevant_IVAR_rows = self.IVAR_rows[relevant_QSOs,first_relevant_cell:last_relevant_cell+1]
        relevant_LOGLAM_MAP = self.LOGLAM_MAP[first_relevant_cell:last_relevant_cell+1]

        #For each cell, determine the number of skewers for which it is relevant.
        N_relevant_skewers = np.sum(relevant_IVAR_rows,axis=0)
        relevant_cells = N_relevant_skewers>0

        #Calculate F_DELTA_rows from F_rows.
        #Introduce a small 'hack' in order to get around the problem of having cells with no skewers contributing to them.
        # TODO: find a neater way to deal with this
        small = 1.0e-10
        relevant_F_BAR = np.average(relevant_F_rows,weights=relevant_IVAR_rows+small,axis=0)
        relevant_F_DELTA_rows = ((relevant_F_rows)/relevant_F_BAR - 1)*relevant_IVAR_rows

        #Calculate the mean in each cell of the gaussian delta and its square.
        GDB = np.average(relevant_GAUSSIAN_DELTA_rows,weights=relevant_IVAR_rows+small,axis=0)*relevant_cells
        GDSB = np.average(relevant_GAUSSIAN_DELTA_rows**2,weights=relevant_IVAR_rows+small,axis=0)*relevant_cells

        #Calculate the mean in each cell of the density delta and its square.
        DDB = np.average(relevant_DENSITY_DELTA_rows,weights=relevant_IVAR_rows+small,axis=0)*relevant_cells
        DDSB = np.average(relevant_DENSITY_DELTA_rows**2,weights=relevant_IVAR_rows+small,axis=0)*relevant_cells

        #Calculate the mean in each cell of the flux and its square.
        FB = np.average(relevant_F_rows,weights=relevant_IVAR_rows+small,axis=0)*relevant_cells
        FSB = np.average(relevant_F_rows**2,weights=relevant_IVAR_rows+small,axis=0)*relevant_cells

        #Calculate the mean in each cell of the flux delta and its square.
        FDB = np.average(relevant_F_DELTA_rows,weights=relevant_IVAR_rows+small,axis=0)*relevant_cells
        FDSB = np.average(relevant_F_DELTA_rows**2,weights=relevant_IVAR_rows+small,axis=0)*relevant_cells

        #Stitch together the means into a binary table.
        dtype = [('N', '>f4'),('GAUSSIAN_DELTA', '>f4'), ('GAUSSIAN_DELTA_SQUARED', '>f4'), ('DENSITY_DELTA', '>f4'), ('DENSITY_DELTA_SQUARED', '>f4')
                , ('F', '>f4'), ('F_SQUARED', '>f4'), ('F_DELTA', '>f4'), ('F_DELTA_SQUARED', '>f4')]
        statistics = np.array(list(zip(N_relevant_skewers,GDB,GDSB,DDB,DDSB,FB,FSB,FDB,FDSB)),dtype=dtype)

        return statistics


    """
    THE FUNCTIONS BELOW THIS POINT ARE CURRENTLY UNUSED, AND ARE NOT EXPECTED TO BE USED IN FUTURE.

    # TODO: get rid of this function as the one below now covers it
    #Method to extract all data from an input file of a given format.
    @classmethod
    def get_all_data(cls,filename,file_number,input_format,lambda_min=0,IVAR_cutoff=lya,SIGMA_G=None,gaussian_only=False):

        lya = 1215.67
        h = fits.open(filename)

        h_R, h_Z, h_D, h_V = get_COSMO(h,input_format)
        h_lya_lambdas = get_lya_lambdas(h,input_format)

        #Calculate the first_relevant_cell.
        first_relevant_cell = np.argmax(h_lya_lambdas >= lambda_min)

        if input_format == 'physical_colore':

            #Extract data from the HDUlist.
            TYPE = h[1].data['TYPE']
            RA = h[1].data['RA']
            DEC = h[1].data['DEC']
            Z_QSO = h[1].data['Z_COSMO']
            DZ_RSD = h[1].data['DZ_RSD']
            DENSITY_DELTA_rows = h[2].data[:,first_relevant_cell:]
            VEL_rows = h[3].data[:,first_relevant_cell:]
            Z = h[4].data['Z'][first_relevant_cell:]
            R = h[4].data['R'][first_relevant_cell:]
            D = h[4].data['D'][first_relevant_cell:]
            V = h[4].data['V'][first_relevant_cell:]

            #Derive the number of quasars and cells in the file.
            N_qso = RA.shape[0]
            N_cells = Z.shape[0]
            if SIGMA_G == None:
                SIGMA_G = h[4].header['SIGMA_G']

            #Derive the MOCKID and LOGLAM_MAP.
            MOCKID = get_MOCKID(h,input_format,file_number)
            LOGLAM_MAP = np.log10(lya*(1+Z))

            #Calculate the Gaussian skewers.
            GAUSSIAN_DELTA_rows = lognormal_to_gaussian(DENSITY_DELTA_rows,SIGMA_G,D)

            if not gaussian_only:
                #Calculate the transmitted flux.
                A,ALPHA,TAU_rows = get_tau(Z,DENSITY_DELTA_rows+1)
                F_rows = np.exp(-TAU_rows)
            else:
                A = None
                ALPHA = None
                TAU_rows = None
                F_rows = None
                DENSITY_DELTA_rows = None

            #Insert placeholder values for remaining variables.
            PLATE = MOCKID
            MJD = np.zeros(N_qso)
            FIBER = np.zeros(N_qso)

            IVAR_rows = make_IVAR_rows(IVAR_cutoff,Z_QSO,LOGLAM_MAP)

            # TODO: Think about how to do this. Also make sure to implement everwhere!
            #Construct grouping variables for appearance.
            #I =
            #II =
            #III =
            #IV =

        elif input_format == 'gaussian_colore':

            #Extract data from the HDUlist.
            TYPE = h[1].data['TYPE']
            RA = h[1].data['RA']
            DEC = h[1].data['DEC']
            Z_QSO = h[1].data['Z_COSMO']
            DZ_RSD = h[1].data['DZ_RSD']
            GAUSSIAN_DELTA_rows = h[2].data[:,first_relevant_cell:]
            VEL_rows = h[3].data[:,first_relevant_cell:]
            Z = h[4].data['Z'][first_relevant_cell:]
            R = h[4].data['R'][first_relevant_cell:]
            D = h[4].data['D'][first_relevant_cell:]
            V = h[4].data['V'][first_relevant_cell:]

            #Derive the number of quasars and cells in the file.
            N_qso = RA.shape[0]
            N_cells = Z.shape[0]
            if SIGMA_G == None:
                SIGMA_G = h[4].header['SIGMA_G']

            #Derive the MOCKID and LOGLAM_MAP.
            MOCKID = get_MOCKID(h,input_format,file_number)
            LOGLAM_MAP = np.log10(lya*(1+Z))

            if not gaussian_only:
                #Calculate the Gaussian skewers.
                DENSITY_DELTA_rows = gaussian_to_lognormal(GAUSSIAN_DELTA_rows,SIGMA_G,D)

                #Calculate the transmitted flux.
                A,ALPHA,TAU_rows = get_tau(Z,DENSITY_DELTA_rows+1)
                F_rows = np.exp(-TAU_rows)
            else:
                A = None
                ALPHA = None
                TAU_rows = None
                F_rows = None
                DENSITY_DELTA_rows = None

            #Insert placeholder values for remaining variables.
            PLATE = MOCKID
            MJD = np.zeros(N_qso)
            FIBER = np.zeros(N_qso)

            IVAR_rows = make_IVAR_rows(IVAR_cutoff,Z_QSO,LOGLAM_MAP)

        elif input_format == 'picca':

            #Extract data from the HDUlist.
            DENSITY_DELTA_rows = h[0].data.T[:,first_relevant_cell:]
            IVAR_rows = h[1].data.T[:,first_relevant_cell:]
            LOGLAM_MAP = h[2].data[first_relevant_cell:]
            RA = h[3].data['RA']
            DEC = h[3].data['DEC']
            Z_QSO = h[3].data['Z']
            PLATE = h[3].data['PLATE']
            MJD = h[3].data['MJD']
            FIBER = h[3].data['FIBER']
            MOCKID = h[3].data['THING_ID']

            #Derive the number of quasars and cells in the file.
            N_qso = RA.shape[0]
            N_cells = LOGLAM_MAP.shape[0]
            if SIGMA_G == None:
                SIGMA_G = h[4].header['SIGMA_G']

            #Derive Z.
            Z = (10**LOGLAM_MAP)/lya - 1

            #Calculate the Gaussian skewers.
            GAUSSIAN_DELTA_rows = lognormal_to_gaussian(DENSITY_DELTA_rows,SIGMA_G,D)

            if not gaussian_only:
                #Calculate the transmitted flux.
                A,ALPHA,TAU_rows = get_tau(Z,DENSITY_DELTA_rows+1)
                F_rows = np.exp(-TAU_rows)
            else:
                A = None
                ALPHA = None
                TAU_rows = None
                F_rows = None
                DENSITY_DELTA_rows = None

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

        return cls(N_qso,N_cells,SIGMA_G,ALPHA,TYPE,RA,DEC,Z_QSO,DZ_RSD,MOCKID,PLATE,MJD,FIBER,GAUSSIAN_DELTA_rows,DENSITY_DELTA_rows,VEL_rows,IVAR_rows,F_rows,R,Z,D,V,LOGLAM_MAP,A)



    #Method to create a new object from an existing one, having specified which MOCKIDs we want to include.
    # TODO: add something to check that we can just take values from 1 of the objects
    @classmethod
    def choose_qsos(cls,object_A,MOCKIDs):

        rows = ['']*len(MOCKIDs)
        s = set(MOCKIDs)
        j = 0
        for i, qso in enumerate(object_A.MOCKID):
            if qso in s:
                rows[j] = i
                j=j+1

        N_qso = len(rows)
        N_cells = object_A.N_cells
        SIGMA_G = object_A.SIGMA_G
        ALPHA = object_A.ALPHA

        TYPE = object_A.TYPE[rows]
        RA = object_A.RA[rows]
        DEC = object_A.DEC[rows]
        Z_QSO = object_A.Z_QSO[rows]
        DZ_RSD = object_A.DZ_RSD[rows]
        MOCKID = object_A.MOCKID[rows]
        PLATE = object_A.PLATE[rows]
        MJD = object_A.MJD[rows]
        FIBER = object_A.FIBER[rows]

        GAUSSIAN_DELTA_rows = object_A.GAUSSIAN_DELTA_rows[rows,:]
        DENSITY_DELTA_rows = object_A.DENSITY_DELTA_rows[rows,:]
        VEL_rows = object_A.VEL_rows[rows,:]
        IVAR_rows = object_A.IVAR_rows[rows,:]
        F_rows = object_A.F_rows[rows,:]

        Z = object_A.Z
        LOGLAM_MAP = object_A.LOGLAM_MAP
        R = object_A.R
        D = object_A.D
        V = object_A.V
        A = object_A.A

        return cls(N_qso,N_cells,SIGMA_G,ALPHA,TYPE,RA,DEC,Z_QSO,DZ_RSD,MOCKID,PLATE,MJD,FIBER,GAUSSIAN_DELTA_rows,DENSITY_DELTA_rows,VEL_rows,IVAR_rows,F_rows,R,Z,D,V,LOGLAM_MAP,A)

    #Method to create a new object from an existing one, having specified which cells we want to include.
    # TODO: change this so you can specify a z_min/max, or lambda_min/max, rather than just any list of cells. Would need to deal with the case of both z and lambda limits being set.
    @classmethod
    def choose_cells(cls,object_A,cells):

        N_qso = object_A.N_qso
        N_cells = len(cells)
        SIGMA_G = object_A.SIGMA_G
        ALPHA = object_A.ALPHA

        TYPE = object_A.TYPE
        RA = object_A.RA
        DEC = object_A.DEC
        Z_QSO = object_A.Z_QSO
        DZ_RSD = object_A.DZ_RSD
        MOCKID = object_A.MOCKID
        PLATE = object_A.PLATE
        MJD = object_A.MJD
        FIBER = object_A.FIBER

        GAUSSIAN_DELTA_rows = object_A.GAUSSIAN_DELTA_rows[:,cells]
        DENSITY_DELTA_rows = object_A.DENSITY_DELTA_rows[:,cells]
        VEL_rows = object_A.VEL_rows[:,cells]
        IVAR_rows = object_A.IVAR_rows[:,cells]

        Z = object_A.Z[cells]
        LOGLAM_MAP = object_A.LOGLAM_MAP[cells]
        R = object_A.R[cells]
        D = object_A.D[cells]
        V = object_A.V[cells]
        A = object_A.A[cells]

        return cls(N_qso,N_cells,SIGMA_G,TYPE,RA,DEC,Z_QSO,DZ_RSD,MOCKID,PLATE,MJD,FIBER,GAUSSIAN_DELTA_rows,DENSITY_DELTA_rows,VEL_rows,IVAR_rows,R,Z,D,V,LOGLAM_MAP)


    def save(self,filename,header,output_format):

        success = 0
        while success == 0:
            if output_format == 'colore':

                #Organise the data into colore-format arrays.
                colore_1_data = []
                for i in range(self.N_qso):
                    colore_1_data += [(self.TYPE[i],self.RA[i],self.DEC[i],self.Z_QSO[i],self.DZ_RSD[i],self.MOCKID[i])]

                dtype = [('TYPE', '>f4'), ('RA', '>f4'), ('DEC', '>f4'), ('Z_COSMO', '>f4'), ('DZ_RSD', '>f4'), ('MOCKID', int)]
                colore_1 = np.array(colore_1_data,dtype=dtype)

                colore_2 = self.DENSITY_DELTA_rows
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
                picca_0 = self.DENSITY_DELTA_rows.T
                picca_1 = self.IVAR_rows.T
                picca_2 = self.LOGLAM_MAP

                picca_3_data = []
                for i in range(self.N_qso):
                    picca_3_data += [(self.RA[i],self.DEC[i],self.Z_QSO[i],self.PLATE[i],self.MJD[i],self.FIBER[i],self.MOCKID[i])]

                dtype = [('RA', '>f4'), ('DEC', '>f4'), ('Z', '>f4'), ('PLATE', int), ('MJD', '>f4'), ('FIBER', int), ('MOCKID', int)]
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
    def crop(cls,object_A,MOCKID,cells):

        rows = ['']*len(MOCKID)
        s = set(MOCKID)
        j = 0
        for i, qso in enumerate(object_A.MOCKID):
            if qso in s:
                rows[j] = i
                j=j+1

        N_qso = len(rows)
        N_cells = len(cells)

        TYPE = object_A.TYPE[rows]
        RA = object_A.RA[rows]
        DEC = object_A.DEC[rows]
        Z_QSO = object_A.Z_QSO[rows]
        DZ_RSD = object_A.DZ_RSD[rows]
        MOCKID = object_A.MOCKID[rows]
        PLATE = object_A.PLATE[rows]
        MJD = object_A.MJD[rows]
        FIBER = object_A.FIBER[rows]

        DENSITY_DELTA_rows = object_A.DENSITY_DELTA_rows[rows,:]
        DENSITY_DELTA_rows = DENSITY_DELTA_rows[:,cells]

        VEL_rows = object_A.VEL_rows[rows,:]
        VEL_rows = VEL_rows[:,cells]

        IVAR_rows = object_A.IVAR_rows[rows,:]
        IVAR_rows = IVAR_rows[:,cells]

        Z = object_A.Z[cells]
        LOGLAM_MAP = object_A.LOGLAM_MAP[cells]
        R = object_A.R[cells]
        D = object_A.D[cells]
        V = object_A.V[cells]

        return cls(N_qso,N_cells,TYPE,RA,DEC,Z_QSO,DZ_RSD,MOCKID,PLATE,MJD,FIBER,DENSITY_DELTA_rows,VEL_rows,IVAR_rows,R,Z,D,V,LOGLAM_MAP)

    #NOT YET READY TO BE USED - maybe not necessary?
    #Method to extract reduced data from a set of input files of a given format, with a given list of MOCKIDs for each file.
    @classmethod
    def WIP(cls,file_infos,input_format,z_min):


        for file_info in file_infos:

            lya = 1215.67

            h = fits.open(file_info[filename])

            h_MOCKID = get_MOCKID(h,input_format,file_info[file_number])
            h_R, h_Z, h_D, h_V = get_COSMO(h,input_format)

            #Work out which rows in the hdulist we are interested in.
            rows = ['']*len(file_info[MOCKIDs])
            s = set(file_info[MOCKIDs])
            j = 0
            for i, qso in enumerate(h_MOCKID):
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
            DENSITY_DELTA_rows = []
            VEL_rows = []
            Z = []
            R = []
            D = []
            V = []
            N_qso = []
            N_cells = []
            MOCKID = []
            LOGLAM_MAP = []
            PLATE = []
            MJD = []
            FIBER = []
            IVAR_rows = []

            if input_format == 'physical_colore':

                #Extract data from the HDUlist.
                TYPE = np.concatenate((TYPE,h[1].data['TYPE'][rows]))
                RA = np.concatenate((RA,h[1].data['RA'][rows]))
                DEC = np.concatenate((DEC,h[1].data['DEC'][rows]))
                Z_QSO = np.concatenate((Z_QSO,h[1].data['Z_COSMO'][rows]))
                DZ_RSD = np.concatenate((DZ_RSD,h[1].data['DZ_RSD'][rows]))

                DENSITY_DELTA_rows = np.concatenate((DENSITY_DELTA_rows,h[2].data[rows,first_relevant_cell:]),axis=0)

                VEL_rows = np.concatenate((VEL_rows,h[3].data[rows,first_relevant_cell:]),axis=0)

                Z = h[4].data['Z'][first_relevant_cell:]
                R = h[4].data['R'][first_relevant_cell:]
                D = h[4].data['D'][first_relevant_cell:]
                V = h[4].data['V'][first_relevant_cell:]

                #Derive the number of quasars and cells in the file.
                N_qso = N_qso + RA.shape[0]
                N_cells = Z.shape[0]

                #Derive the MOCKID and LOGLAM_MAP.
                MOCKID = np.concatenate((MOCKID,h_MOCKID[rows]))
                LOGLAM_MAP = np.log10(lya*(1+Z))

                #Insert placeholder values for remaining variables.
                PLATE = np.concatenate((PLATE,np.zeros(N_qso)))
                MJD = np.concatenate((MJD,np.zeros(N_qso)))
                FIBER = np.concatenate((FIBER,np.zeros(N_qso)))
                IVAR_rows = np.concatenate((IVAR_rows,np.ones((N_qso,N_cells))),axis=0)

            elif input_format == 'picca':

                #Extract data from the HDUlist.
                DENSITY_DELTA_rows = h[0].data.T[rows,first_relevant_cell:]

                IVAR_rows = h[1].data.T[rows,first_relevant_cell:]

                LOGLAM_MAP = h[2].data[first_relevant_cell:]

                RA = h[3].data['RA'][rows]
                DEC = h[3].data['DEC'][rows]
                Z_QSO = h[3].data['Z'][rows]
                PLATE = h[3].data['PLATE'][rows]
                MJD = h[3].data['MJD'][rows]
                FIBER = h[3].data['FIBER'][rows]
                MOCKID = h[3].data['MOCKID'][rows]

                #Derive the number of quasars and cells in the file.
                N_qso = RA.shape[0]
                N_cells = LOGLAM_MAP.shape[0]

                #Derive Z.
                Z = (10**LOGLAM_MAP)/lya - 1

                #Insert placeholder variables for remaining variables.
                TYPE = np.zeros(RA.shape[0])
                R = np.zeros(Z.shape[0])
                D = np.zeros(Z.shape[0])
                V = np.zeros(Z.shape[0])
                DZ_RSD = np.zeros(RA.shape[0])
                VEL_rows = np.zeros(DENSITY_DELTA_rows.shape)

            else:
                print('Input format not recognised: current options are "colore" and "picca".')
                print('Please choose one of these options and try again.')

            h.close()

        return cls(N_qso,N_cells,TYPE,RA,DEC,Z_QSO,DZ_RSD,MOCKID,PLATE,MJD,FIBER,DENSITY_DELTA_rows,VEL_rows,IVAR_rows,R,Z,D,V,LOGLAM_MAP)
    """
