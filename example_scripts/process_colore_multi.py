import numpy as np
from astropy.io import fits
import process_functions as functions
from multiprocessing import Pool
import sys
import time

start = time.time()

#Top level script to manage the conversion of CoLoRe output files into files more useful for analysis.
#These are produced on a per-HEALPix pixel basis.

lya = 1215.67

#Define the desired power of 2 for Nside for the output. This should be larger than that of the input.
N_side_pow2 = 3
N_side = 2**N_side_pow2
N_pix = 12*N_side**2

#Define the original file structure
original_file_location = '/Users/jfarr/Projects/repixelise/test_input/'
original_filename_structure = 'N10000_out_srcs_s1_{}.fits' #file_number
file_numbers = [16]
input_format = 'physical_colore'

#Set file structure
new_base_file_location = '/Users/jfarr/Projects/repixelise/test_output/test_multi/'
new_file_structure = '{}/{}/'               #pixel number//100, pixel number
new_filename_structure = '{}-{}-{}.fits'    #file type, nside, pixel number

#Choose options
lambda_min = 3550
enforce_picca_zero_mean_delta = True

#Calculate the minimum value of z that we are interested in.
#i.e. the z value for which lambda_min cooresponds to the lya wavelength.
z_min = lambda_min/lya - 1

#Get the simulation parameters from the parameter file.
parameter_filename = 'out_params.cfg'
simulation_parameters = functions.get_simulation_parameters(original_file_location,parameter_filename)

# TODO: Modify this to accomodate other density types.
#If density type is not lognormal, then crash.
if simulation_parameters['dens_type'] != 0:
    error('Density is not lognormal. Non-lognormal densities are not currently supported.')



#Make master file
master_data, bad_coordinates_data, file_pixel_map = functions.get_ID_data(original_file_location,original_filename_structure,input_format,file_numbers,N_side)

#Write master file.
master_filename = new_base_file_location + '/' + 'nside_{}_'.format(N_side) + 'master.fits'
functions.write_ID(master_filename,master_data,N_side)
print('Master file containing {} objects saved at :\n{}\n'.format(master_data.shape[0],master_filename))

#Write bad coordinates file.
bad_coordinates_filename = new_base_file_location + '/' + 'nside_{}_'.format(N_side) + 'bad_coordinates.fits'
functions.write_ID(bad_coordinates_filename,bad_coordinates_data,N_side)
print('"Bad coordinates" file containing {} objects saved at:\n{}\n'.format(bad_coordinates_data.shape[0],bad_coordinates_filename))



#Make a list of the pixels that the files cover.
pixel_list = list(sorted(set(master_data['PIXNUM'])))
file_number_list = list(sorted(set([int(functions.number_to_string(qso['MOCKID'],10)[:-7]) for qso in master_data])))


print('Making master file: {}s.'.format(time.time()-start))
start=time.time()


#Make the new file structure
functions.make_file_structure(new_base_file_location,pixel_list)

#Define the pixelisation process.
#For each pixel, create a 'pixel object' by using the master file's data to collect QSOs from the same pixel, stored in different files.
#Then save this object in the various formats required.
print(' \n')
def pixelise(pixel,original_file_location,original_filename_structure,input_format,master_data,pixel_list,file_number_list,file_pixel_map,z_min,new_base_file_location,new_file_structure,N_side):

    print('Working on pixel number {}.'.format(pixel))
    location = new_base_file_location + new_file_structure.format(pixel//100,pixel)

    #Make file into an object
    pixel_object = functions.make_pixel_object(pixel,original_file_location,original_filename_structure,input_format,master_data,pixel_list,file_number_list,file_pixel_map,z_min)

    #Make some useful headers
    header = fits.Header()
    header['NSIDE'] = N_side
    header['PIXNUM'] = pixel
    header['LYA'] = lya
    header['NQSO'] = pixel_object.N_qso

    #Gaussian CoLoRe
    filename = new_filename_structure.format('gaussian-colore',N_side,pixel)
    pixel_object.save_as_gaussian_colore(location,filename,header)

    #lognorm CoLoRe
    filename = new_filename_structure.format('physical-colore',N_side,pixel)
    pixel_object.save_as_physical_colore(location,filename,header)

    #Picca density
    filename = new_filename_structure.format('picca-density',N_side,pixel)
    pixel_object.save_as_picca_density(location,filename,header,zero_mean_delta=True*enforce_picca_zero_mean_delta,lambda_min=lambda_min)

    #transmission
    filename = new_filename_structure.format('transmission',N_side,pixel)
    pixel_object.save_as_transmission(location,filename,header,lambda_min=lambda_min)

    #picca flux
    filename = new_filename_structure.format('picca-flux',N_side,pixel)
    pixel_object.save_as_picca_flux(location,filename,header,zero_mean_delta=True*enforce_picca_zero_mean_delta,lambda_min=lambda_min)

    #Close everything
    print('\nPixel number {} complete!\n'.format(pixel))

    return

#Use multiprocessing to pixelise quickly.
N_processes = int(sys.argv[1])
pool = Pool(processes = N_processes)

tasks = []
for pixel in pixel_list:
    tasks += [(pixel,original_file_location,original_filename_structure,input_format,master_data,pixel_list,file_number_list,file_pixel_map,z_min,new_base_file_location,new_file_structure,N_side)]

if __name__ == '__main__':
    pool = Pool(processes = N_processes)
    results = pool.starmap  (pixelise,tasks)

print('Making pixel files: {}s.'.format(time.time()-start))
