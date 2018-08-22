import numpy as np
from astropy.io import fits
import time
import os
import healpy as hp

lya_rest = 1215.67

#Define a function to print a progress bar.
def progress_bar(N_complete,N_tasks,start_time):

    N_chunks = 20
    N_chunks_complete = int((N_complete*N_chunks) // (N_tasks))
    block_char = '-'
    progress_bar = '|' + block_char*N_chunks_complete + ' '*(N_chunks-N_chunks_complete) + '|'

    current_time = time.time()
    time_elapsed = current_time - start_time
    estimated_time_remaining = (time_elapsed)*((N_tasks-N_complete)/N_complete)
    print(' -> current progress: {} {:4d} of {:4d} complete ({:3.0%}), {:4.0f}s elapsed, ~{:5.0f}s remaining'.format(progress_bar,N_complete,N_tasks,N_complete/N_tasks,time_elapsed,estimated_time_remaining),end="\r")

    if N_complete == N_tasks:
        print('\nProcess complete!')

    return

#Function to create a file structure based on a set of numbers, of the form "x//100"/"x".
def make_file_structure(base_location,numbers):

    first_level = []
    for number in numbers:
        first_level += [number//100]

    first_level_set = list(sorted(set(first_level)))

    for i in first_level_set:

        os.mkdir(base_location+'/'+str(i))

        for j, number in enumerate(numbers):

            if first_level[j] == i:
                os.mkdir(base_location+'/'+str(i)+'/'+str(number))

    return

# NOTE: don't think this is used
#Function to convert a list of numbers to a list of n-digit strings.
def number_to_string(number,string_length):

    number_str = str(number)
    if len(number_str)<=string_length:
        string = '0'*(string_length-len(number_str))+number_str
    else:
        exit('The file number is too great to construct a unique MOCKID (more than 3 digits).')

    return string

#Function to construct an array of MOCKIDs given a file_number and a list of row_numbers.
def make_MOCKID(file_number,row_numbers):

    N_qso = len(row_numbers)
    node = '0'*(len(str(file_number))-5) + str(file_number)

    MOCKID = ['']*N_qso
    for i in range(N_qso):
        row_numbers[i] = str(row_numbers[i]+1)
        if len(row_numbers[i])<=7:
            row_numbers[i] = '0'*(7-len(row_numbers[i]))+row_numbers[i]
        else:
            exit('The row number is too great to construct a unique MOCKID (more than 7 digits).')
        MOCKID[i] = int(node+row_numbers[i])

    MOCKID = np.array(MOCKID)

    return MOCKID

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
            pixel_ID[i] = int(hp.pixelfunc.ang2pix(N_side,theta[i],phi[i],nest=True))
        else:
            pixel_ID[i] = -1

    return pixel_ID

#Function to make ivar mask
def make_IVAR_rows(IVAR_cutoff,Z_QSO,LOGLAM_MAP):

    N_cells = LOGLAM_MAP.shape[0]
    N_qso = Z_QSO.shape[0]

    lya_lambdas = IVAR_cutoff*(1+Z_QSO)
    IVAR_rows = np.ones((N_qso,N_cells),dtype='float32')
    lambdas = 10**LOGLAM_MAP

    for i in range(N_qso):
        last_relevant_cell = np.searchsorted(lambdas,lya_lambdas[i]) - 1

        for j in range(last_relevant_cell+1,N_cells):
            IVAR_rows[i,j] = 0.

    return IVAR_rows

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

#Function to retrieve relevant simulation parameters from the param.cfg file.
def get_simulation_parameters(location,filename):

    #Create a string of the parameter file to search.
    parameter_str = open(location + '/' + filename,'r').read()

    #Define the parameters to search for and the intricacies of the parameter file.
    divider = '\n'
    #parameters = [('dens_type','int'),('r_smooth','f4'),('n_grid','int'),('gaussian_skewers','str'),('omega_M','f4'),('omega_L','f4'),('omega_B','f4'),('h','f4'),('w','f4'),('ns','f4'),('sigma_8','f4')]
    parameters = [('dens_type','int'),('r_smooth','f4'),('n_grid','int'),('omega_M','f4'),('omega_L','f4'),('omega_B','f4'),('h','f4'),('w','f4'),('ns','f4'),('sigma_8','f4')]
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
def normalise_deltas(DELTA_rows,mean_DELTA):

    N_cells = DELTA_rows.shape[1]
    DELTA_rows_normalised = np.zeros(DELTA_rows.shape)

    for j in range(N_cells):
        DELTA_rows_normalised[:,j] = (DELTA_rows[:,j] + 1)/(mean_DELTA[j] + 1) - 1

    return DELTA_rows_normalised

#Function to interpolate via the NGP method.
def get_NGPs(x,x_new):

    NGPs = np.zeros(x_new.shape)

    for i,x_new_value in enumerate(x_new):
        distances2 = (x-x_new_value)**2
        NGPs[i] = np.argmin(distances2)

    return NGPs

#Function to return the index of the point in a sorted array closest to a given value.
def NN_sorted(arr,val):

    N = arr.shape[0]
    i = np.searchsorted(arr,val)

    if i > 0 and i < N:
        if abs(arr[i] - val) > abs(arr[i-1] - val):
            i -= 1
    elif i == N:
        i -= 1

    return i

#transform differential length (in Mpc/h) to differential velocity (in km/s)
def get_dkms_dhMpc(z,Om=0.3147):

    E_z = np.sqrt(Om*(1+z)**3 + (1-Om))
    dkms_dhMpc = 100. * E_z / (1+z)

    return dkms_dhMpc

#Function to check if two objects are the same, return an error if not.
def confirm_identical(A,B,item_name=None,array=False):
    if not array:
        if A != B:
            if item_name:
                raise ValueError('Wrongly attempted to combine non-identical objects of type {}!'.format(item_name))
            else:
                raise ValueError('Wrongly attempted to combine non-identical objects!')
            result = False
        else:
            result = True
    else:
        if any(A != B):
            if item_name:
                raise ValueError('Wrongly attempted to combine non-identical objects of type {}!'.format(item_name))
            else:
                raise ValueError('Wrongly attempted to combine non-identical objects!')
            result = False
        else:
            result = True
    return result
