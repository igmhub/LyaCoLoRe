import numpy as np
from astropy.io import fits
import time
import os
import healpy as hp
import shutil
import gzip
import sys
import warnings

lya_rest = 1215.67

#Function to get the file name using a defined structure.
def get_file_name(base_dir,base_name,nside,pixel,compressed=False):
    if compressed:
        return base_dir+'/{}-{}-{}.fits.gz'.format(base_name,nside,pixel)
    else:
        return base_dir+'/{}-{}-{}.fits'.format(base_name,nside,pixel)

#Function to get the directory name using a defined structure.
def get_dir_name(base_dir,pixel):
    return base_dir+'/{}/{}/'.format(pixel//100,pixel)

def compress_file(filename,ext='.gz',remove=True):
    with open(filename,'rb') as f_in, gzip.open(filename+ext,'wb') as f_out:
        shutil.copyfileobj(f_in,f_out)
    if remove:
        os.remove(filename)
    return

def get_edges(x):
    x_edges = np.concatenate([[(x[0]-(x[1]-x[0])/2)],(x[1:]+x[:-1])/2,[(x[-1]+(x[-1]-x[-2])/2)]])
    return x_edges

def get_centres(x_edges):
    x = (x_edges[1:] + x_edges[:-1])/2.
    return x

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

    numbers = np.array(numbers)
    first_level = numbers//100
    first_level_set = list(sorted(set(first_level)))

    for i in first_level_set:
        try:
            os.mkdir(base_location+'/'+str(i))
        except FileExistsError:
            pass

    for n in numbers:
        try:
            os.mkdir(base_location+'/'+str(n//100)+'/'+str(n))
        except FileExistsError:
            pass

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
    #We check that the angular coordinates are valid.
    #Give all objects with invalid coordinates an "error" ID number (-1).
    valid_QSOs = (0 <= theta) * (theta <= np.pi) * (0 <= phi) * (phi <= 2*np.pi)
    pixel_ID = np.ones_like(RA) * (-1)
    pixel_ID[valid_QSOs] = hp.pixelfunc.ang2pix(N_side,theta[valid_QSOs],phi[valid_QSOs],nest=True)

    return pixel_ID

#Function to return the neighbouring HEALPix pixels to a given HEALPix pixel.
def get_pixel_neighbours(pixel,N_side=16):
    theta,phi = hp.pix2ang(N_side,pixel,nest=True)
    neighbours = hp.get_all_neighbours(N_side,theta=theta,phi=phi,nest=True)
    return neighbours

#Function to add neighbouring HEALPix pixels to an array of HEALPix pixels.
def add_pixel_neighbours(pixels,N_side=16):
    final_pixel_set = set(pixels)
    for pixel in pixels:
        neighbours = get_pixel_neighbours(pixel,N_side=N_side)
        final_pixel_set = final_pixel_set.union(set(neighbours))
    return np.array(list(final_pixel_set))

#TODO add pixel_list functionality.
#Function to return a filter function for restricting the QSO footprint.
def make_QSO_filter(footprint,N_side=16,pixel_list=None):

    #See if we can use desimodel. This is preferable as it will be the most
    #up-do-date footprint.
    try:
        from desimodel.footprint import tiles2pix, is_point_in_desi
        from desimodel.io import load_tiles
        desimodel_installed = True
    except ModuleNotFoundError:
        desimodel_installed = False

    #If we have desimodel and want to replicate the footprint precisely, use
    #function "is_point_in_desi".
    if desimodel_installed and footprint=='desi':
        tiles = load_tiles()
        def QSO_filter(RA,DEC):
            return is_point_in_desi(tiles,RA,DEC)

    #If not, but we still want to filter...
    elif footprint in ['desi','desi_pixel','desi_pixel_plus']:

        #If desimodel is installed, then we use "tiles2pix" to determine which
        #pixels to include.
        if desimodel_installed:
            from desimodel.footprint import tiles2pix
            if footprint=='desi_pixel':
                valid_pixels = tiles2pix(N_side)
            elif footprint=='desi_pixel_plus':
                valid_pixels = tiles2pix(N_side)
                valid_pixels = add_pixel_neighbours(valid_pixels)

        #Otherwise, we load pixel lists from file. Note: using desimodel is
        #preferable to loading from file as the footprint could change, and
        #desimodel will be more up to date than the lists in this case.
        else:
            print('desimodel not installed; loading pixel footprints from file...')
            if footprint=='desi':
                print('desi footprint not possible without desimodel: using desi_pixel instead...')
                valid_pixels = np.loadtxt('input_files/pixel_footprints/DESI_pixels.txt',dtype=int)
            elif footprint=='desi_pixel':
                valid_pixels = np.loadtxt('input_files/pixel_footprints/DESI_pixels.txt',dtype=int)
            elif footprint=='desi_pixel_plus':
                valid_pixels = np.loadtxt('input_files/pixel_footprints/DESI_pixels_plus.txt',dtype=int)

        #With a list of valid pixels, we now can make a filter.
        def QSO_filter(RA,DEC):
            theta = (np.pi/180.0)*(90.0-DEC)
            phi = (np.pi/180.0)*RA
            pix = hp.ang2pix(N_side,theta,phi,nest=True)
            w = np.in1d(pix,valid_pixels)
            return w

    #Else if we don't want to filter at all, set the filter to "None".
    elif footprint=='full_sky':
        def QSO_filter(RA,DEC):
            return np.ones(RA.shape).astype('bool')

    else:
        print('Footprint not recognised; no filter applied.')
        def QSO_filter(RA,DEC):
            return np.ones(RA.shape).astype('bool')

    return QSO_filter

#Function to return a filter for restricting the QSO footprint.
def choose_filter(desi_footprint,desi_footprint_pixel,desi_footprint_pixel_plus,desimodel_installed,N_side=16,pixel_list=None):

    if np.sum((desi_footprint,desi_footprint_pixel,desi_footprint_pixel_plus)) > 1:
            raise ValueError('Please choose only 1 type of DESI footprint.')

    if desimodel_installed:
        from desimodel.footprint import tiles2pix, is_point_in_desi
        if desi_footprint:
            def QSO_filter(RA,DEC):
                return is_point_in_desi(tiles,RA,DEC)
        elif desi_footprint_pixel:
            QSO_filter = tiles2pix(N_side)
        elif desi_footprint_pixel_plus:
            QSO_filter = tiles2pix(N_side)
            QSO_filter = add_pixel_neighbours(QSO_filter)
        else:
            QSO_filter = None
    else:
        if desi_footprint or desi_footprint_pixel:
            QSO_filter = np.loadtxt('input_files/DESI_pixels.txt',dtype=int)
        elif desi_footprint_pixel_plus:
            QSO_filter = np.loadtxt('input_files/DESI_pixels_plus.txt',dtype=int)
        else:
            QSO_filter = None

    return QSO_filter

#Function to make ivar mask
def make_IVAR_rows(IVAR_cutoff,Z_QSO,LOGLAM_MAP):

    #Make an array of rest frame lambdas.
    lambdas = 10**LOGLAM_MAP
    lambdas_rf = np.outer(1/(1+Z_QSO),lambdas)

    #Filter according to the cutoff.
    IVAR_rows = (lambdas_rf <= IVAR_cutoff).astype('float32')

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

    NGPs = NGPs.astype(int)

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

#Function to rebin skewers according to a fixed cell size.
def rebin_skewers(skewers,R,cell_size):

    """
    grid_start = R[0]
    rebin_map = np.array((R - R[0])//rebin_size_hMpc,dtype='int')
    rebin_new_delta = np.zeros((N_qso,max(rebin_map)))
    j_hi = -1
    for j in range(max(rebin_map)):
        j_lo = j_hi + 1 #np.argmax(rebin_map==j)
        j_hi = np.searchsorted(rebin_map,j+1) - 1
        rebin_new_delta[:,j] = np.average(new_delta[:,j_lo:j_hi+1],axis=1)
    """

    R_edges = R[1:] - R[:-1]

    R_shift = R - R[0]

    R_new = np.arange(R[0],R[-1],cell_size)
    R_new_edges = R_new[1:] - R_new[:-1]

    N_cells = R.shape[0]

    return rebinned_skewers

#Function to rebin skewers by merging a certain number of cells together.
def merge_cells(rows,N_merge):

    #Work out if it's an array or vector. Trim the rows so we have no half full.
    N_dim = len(rows.shape)
    if N_dim == 1:
        N_cells = rows.shape[0]
        N_cells_new = N_cells//N_merge
        rows = rows[:N_cells_new*N_merge]
        new_shape = (N_cells_new,N_merge)
        merged_rows = rows.reshape(new_shape).mean(-1)
    elif N_dim == 2:
        N_rows = rows.shape[0]
        N_cells = rows.shape[1]
        N_cells_new = N_cells//N_merge
        rows = rows[:,:N_cells_new*N_merge]
        new_shape = (N_rows,N_cells_new,N_merge)
        merged_rows = rows.reshape(new_shape).mean(-1)

    """
    #Fill in the shell with merged values
    for j in range(N_cells_new):
        if N_dim == 1:
            merged_rows[j] = np.average(rows[j*N_merge:(j+1)*N_merge])
        elif N_dim == 2:
            merged_rows[:,j] = np.average(rows[:,j*N_merge:(j+1)*N_merge],axis=1)
    """
    return merged_rows

#Function to renormalise and rebin data in a picca file.
def renorm_rebin_picca_file(filepath,old_mean=None,new_mean=None,N_merge=None,IVAR_cutoff=1150.,min_number_cells=2,out_filepath=None,overwrite=False,compress=False):

    #Open the existing file.
    h = fits.open(filepath)
    old_delta_rows = h[0].data.T
    header = h[0].header
    hdu_iv = h[1]
    hdu_LOGLAM_MAP = h[2]
    hdu_CATALOG = h[3]

    #print('--> renorm-ing...')
    t = time.time()
    if old_mean is None and new_mean is None:
        renorm_delta_rows = old_delta_rows
    elif old_mean is None and new_mean is not None:
        cells = new_mean>0
        renorm_delta_rows = np.zeros(old_delta_rows.shape)
        renorm_delta_rows[:,cells] = old_delta_rows[:,cells]/new_mean[cells] - 1
    elif old_mean is not None and new_mean is None:
        renorm_delta_rows = (1 + old_delta_rows) * old_mean
    else:
        #Renormalise the deltas.
        cells = new_mean>0
        renorm_delta_rows = np.zeros(old_delta_rows.shape)
        renorm_delta_rows[:,cells] = (old_mean[cells]*(1 + old_delta_rows[:,cells]))/new_mean[cells] - 1
    #print('----> {:1.3f}s'.format(time.time()-t))

    #print('--> rebinning...')
    t = time.time()
    #Rebin the deltas if necessary.
    if N_merge is not None:
        if N_merge > 1:
            #Merge the deltas and make a new LOGLAM_MAP.
            new_delta_rows = merge_cells(renorm_delta_rows,N_merge).astype('float32')
            old_lambdas = 10**hdu_LOGLAM_MAP.data
            new_LOGLAM_MAP = np.log10(merge_cells(old_lambdas,N_merge)).astype('float32')

            #Determine which new cells are made of entirely
            Z_QSO = hdu_CATALOG.data['Z']
            new_iv_rows = (merge_cells(hdu_iv.data.T,N_merge)==1).astype('float32')
            #new_iv_rows = make_IVAR_rows(IVAR_cutoff,Z_QSO,new_LOGLAM_MAP)
            relevant_QSOs = (np.sum(new_iv_rows,axis=1)>min_number_cells)
            new_delta_rows = new_delta_rows[relevant_QSOs,:].astype('float32')
            new_iv_rows = new_iv_rows[relevant_QSOs,:].astype('float32')
            catalog_data = hdu_CATALOG.data[relevant_QSOs]

            #Reconstruct the non-delta HDUs.
            hdu_deltas_new = fits.PrimaryHDU(data=new_delta_rows.T,header=header)
            hdu_iv_new = fits.ImageHDU(data=new_iv_rows.T,header=hdu_iv.header,name='IV')
            hdu_LOGLAM_MAP_new = fits.ImageHDU(data=new_LOGLAM_MAP,header=hdu_LOGLAM_MAP.header,name='LOGLAM_MAP')
            hdu_CATALOG_new = fits.BinTableHDU(catalog_data,header=hdu_CATALOG.header,name='CATALOG')

        else:
            new_delta_rows = renorm_delta_rows.astype('float32')
            hdu_deltas_new = fits.PrimaryHDU(data=new_delta_rows.T,header=header)
            hdu_iv_new = hdu_iv
            hdu_LOGLAM_MAP_new = hdu_LOGLAM_MAP
            hdu_CATALOG_new = hdu_CATALOG
    else:
        new_delta_rows = renorm_delta_rows.astype('float32')
        hdu_deltas_new = fits.PrimaryHDU(data=new_delta_rows.T,header=header)
        hdu_iv_new = hdu_iv
        hdu_LOGLAM_MAP_new = hdu_LOGLAM_MAP
        hdu_CATALOG_new = hdu_CATALOG
    #print('----> {:1.3f}s'.format(time.time()-t))

    #print('--> saving...')
    t = time.time()
    #Save the file again.
    hdulist = fits.HDUList([hdu_deltas_new, hdu_iv_new, hdu_LOGLAM_MAP_new, hdu_CATALOG_new])
    if not out_filepath:
        out_filepath = filepath
    hdulist.writeto(out_filepath,overwrite=overwrite)
    hdulist.close()
    h.close()
    #print('----> {:1.3f}s'.format(time.time()-t))

    #print('--> compressing...')
    t = time.time()
    #Compress the file if desired.
    if compress:
        compress_file(out_filepath)
    #print('----> {:1.3f}s'.format(time.time()-t))

    return

#Function to produce values of a quadratic log functional form.
def quadratic_log(x,A0,A1,A2):
    return np.log(A0) + A1*np.log(x) + A2*(np.log(x))**2

#Function to extract file numbers from filenames given the location and format.
def get_file_numbers(original_file_location,input_filename_structure,input_files):

    file_numbers = []

    #We assume that there is some filename structure before and after the file
    #number. We find both before and after parts.
    struc_before_number = input_filename_structure[:input_filename_structure.find('{')]
    struc_after_number = input_filename_structure[input_filename_structure.find('}')+1:]

    for infi in input_files:
        #Strip the location.
        infi = infi[len(original_file_location):]

        #Get rid of any slashes if necessary.
        if '/' in infi:
            infi = infi[len(infi)-infi[::-1].find('/'):]

        #Strip the structure before the file number.
        infi = infi[len(struc_before_number):]

        #Strip the structure after the file number.
        infi = infi[:len(infi)-len(struc_after_number)]

        #Convert to a number and add to the list.
        file_number = int(infi)
        file_numbers += [file_number]

    return file_numbers

def get_zeff(z,rp,rt,weights,rmin=80.,rmax=120.):

    #Find the cells that are within the r range.
    r = np.sqrt(rp**2 + rt**2)
    cells = (r>rmin)*(r<rmax)

    #Average over these cells using the weights.
    zeff = np.average(z[cells],weights=weights[cells])

    return zeff

class Quiet:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        warnings.filterwarnings('ignore', category=UserWarning)

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
        warnings.filterwarnings('default', category=UserWarning)
