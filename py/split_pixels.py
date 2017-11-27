import numpy as np
from astropy.io import fits
from astropy.table import Table
import healpy as hp

#Program to convert a .fits file of HEALPix pixels into several .fits files of smaller HEALPix pixels.
#Takes files of the format outputted by CoLoRe and converts into the same format that picca takes as input.

#Function to take multiple input files, split them each with 'split_file' and create a master.fits file.
def split_pixels(N_side,files,file_numbers,save_location,output_format):

    #Determine number of files to be split.
    N_files = len(files)

    #Run through each file and split it into pixel fits files, adding each pixel's THING_IDs and pixel_IDs to the master.
    for i in range(N_files):
        THING_ID, pixel_ID = split_file(N_side,files[i],file_numbers[i],save_location,output_format)
        THING_ID = list(THING_ID)
        pixel_ID = list(pixel_ID)
        master_ID_data = []

        for j in range(len(THING_ID)):
            master_ID_data += [(THING_ID[j],pixel_ID[j])]


    #Sort the THING_IDs and pixel_IDs into the right order: first by pixel number, and then by THING_ID.
    dtype = [('THING_ID', 'S10'), ('PIX', int)]
    master_ID = np.array(master_ID_data, dtype=dtype)
    master_ID_sort = np.sort(master_ID, order=['PIX','THING_ID'])

    #Take the output from each splitting and make it into a master file
    master_header = fits.Header()
    master_header['NSIDE'] = N_side
    master = fits.BinTableHDU.from_columns(master_ID_sort,header=master_header)

    #Make a primary HDU.
    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)

    #Make the .fits file.
    hdulist = fits.HDUList([prihdu,master])
    hdulist.writeto(save_location + '/' + 'master.fits')

    return


#Function to generate a list of files to be split from the name of the containing directory, the file prefix and the file suffices.
def file_list(file_location,file_prefix,file_suffices):

    files = ['']*len(file_suffices)
    for i in range(len(file_suffices)):
        files[i] = file_location+'/'+file_prefix+str(file_suffices[i])+'.fits'

    return files


#Function to split an individual .fits file into per-pixel .fits files.
#Takes the desried N_side, original .fits file filename and desired save location as inputs.
#Saves the per-pixel .fits files at the specified save location, and outputs THING_ID and pixel_ID lists.
def split_file(N_side,original_filename,file_number,save_location,output_format):
    #Open the original .fits file, and extract the data we need.
    initial = fits.open(original_filename)
    RA = initial[1].data['RA']
    DEC = initial[1].data['DEC']
    z_qso = initial[1].data['Z_COSMO']
    z = initial[4].data['Z']
    DELTA = initial[2].data[:]

    N_cells = z.shape[0]
    N_qso = z_qso.shape[0]
    N_pix = 12*N_side**2

    iv = np.ones((N_qso,N_cells))
    PLATE = np.zeros(N_qso,dtype=np.int)
    MJD = np.zeros(N_qso,dtype=np.int)
    FIBER = np.zeros(N_qso,dtype=np.int)

    #Set THING_ID as a 10 digit string, of which the first 3 digits correspond to the node number, and the last 7 correspond to the row number in the original file.
    THING_ID = ['']*N_qso
    node = str(file_number)
    if len(node)<=3:
        node = '0'*(3-len(node))+node
    else:
        exit('The node number is too great to construct a unique THING_ID (more than 3 digits).')

    row_numbers = list(range(initial[2].data[:].shape[0]))
    for i in range(len(row_numbers)):
        row_numbers[i] = str(row_numbers[i])
        if len(row_numbers[i])<=7:
            row_numbers[i] = '0'*(7-len(row_numbers[i]))+row_numbers[i]
        else:
            exit('The row number is too great to construct a unique THING_ID (more than 7 digits).')
        THING_ID[i] = node+row_numbers[i]

    #Set up the LOGLAM_MAP
    LOGLAM_MAP = np.log10(1215.67*(1+z))

    #Convert the coordinates into new pixel identifier numbers, according to the N_side specified.
    pixel_ID=np.zeros([1,len(RA)])

    #Convert DEC and RA in degrees to theta and phi in radians.
    theta = (np.pi/180.0)*(90.0-DEC)
    phi = (np.pi/180.0)*RA

    #Make a list of the HEALPix pixel coordinate of each quasar.
    for i in range(len(RA)):
        #Can we just import healpy on cori? May need to install
        #Check that the angular coordinates are valid. Put all objects with invalid coordinates into a non-realistic ID number (-1).
        if 0 <= theta[i] <= np.pi and 0 <= phi[i] <= 2*np.pi:
            pixel_ID[0,i] = hp.pixelfunc.ang2pix(N_side,theta[i],phi[i])
        else:
            pixel_ID[0,i] = -int(node)

    #print('There are %d objects with invalid angular coordinates.' % (sum(pixel_ID[i] == -1 for i in range(len(pixel_ID)))))
    print('Details of these objects are stored in a file corresponding to a pixel number of -(node number).')

    #Set up a pixel_ID list to map between objects and their pixel.
    pixel_ID = pixel_ID.reshape(-1)

    #Set up a list of pixels represented in the original .fits file.
    pixel_list = list(np.sort(list(set(pixel_ID))))

    #For each pixel represented in the original .fits file, make a new file.
    for n in pixel_list:

        #Progress check aide.
        if n+1>0:
            print('Working on pixel %d (N_pix = %d)' % (n,N_pix))
            print('There are %d quasars in pixel %d.' % (sum(pixel_ID[i] == n for i in range(len(pixel_ID))),n))
        else:
            print('Working on set of objects with invalid coordinates.')
            print('There are %d quasars with invalid coordinates.' % (sum(pixel_ID[i] == n for i in range(len(pixel_ID)))))

        pixel_indices = [i for i in range(len(pixel_ID)) if pixel_ID[i]==n]

        pixel_DELTA = np.array([DELTA[i,:] for i in pixel_indices])
        pixel_iv = np.array([iv[i,:] for i in pixel_indices])
        pixel_RA = [RA[i] for i in pixel_indices]
        pixel_DEC = [DEC[i] for i in pixel_indices]
        pixel_z_qso = [z_qso[i] for i in pixel_indices]
        pixel_PLATE = [PLATE[i] for i in pixel_indices]
        pixel_MJD = [MJD[i] for i in pixel_indices]
        pixel_FIBER = [FIBER[i] for i in pixel_indices]
        pixel_THING_ID = [THING_ID[i] for i in pixel_indices]

        #Transpose pixel_DELTA and pixel_iv to match picca input.
        pixel_DELTA = np.transpose(pixel_DELTA)
        pixel_iv = np.transpose(pixel_iv)

        if n+1>0:
            print('There are %d quasars in pixel %d.' % (len(pixel_THING_ID),n))
        else:
            print('There are %d quasars with invalid coordinates.' % (len(pixel_THING_ID)))

        if len(pixel_THING_ID)!=0:

            if output_format==0:
                #Make output for from fitsio
                print('Not written yet!')

            if output_format==1:
                #Construct a table for the final hdu.
                col_RA = fits.Column(name='RA', array=pixel_RA, format='E')
                col_DEC = fits.Column(name='DEC', array=pixel_DEC, format='E')
                col_z_qso = fits.Column(name='Z', array=pixel_z_qso, format='E')
                col_PLATE = fits.Column(name='PLATE', array=pixel_PLATE, format='E')
                col_MJD = fits.Column(name='MJD', array=pixel_MJD, format='E')
                col_FIBER = fits.Column(name='FIBER', array=pixel_FIBER, format='E')
                col_THING_ID = fits.Column(name='THING_ID', array=pixel_THING_ID, format='10A')

                cols = fits.ColDefs([col_RA, col_DEC, col_z_qso, col_PLATE, col_MJD, col_FIBER, col_THING_ID])

                #Add a couple of headers to the file.
                header = fits.Header()
                header['NSIDE'] = N_side
                header['NQSO'] = len(pixel_THING_ID)
                header['PIX'] = int(n)

                #Create hdus from the data arrays
                hdu_DELTA = fits.PrimaryHDU(data=pixel_DELTA,header=header)
                hdu_iv = fits.ImageHDU(data=pixel_iv,header=header,name='IV')
                hdu_LOGLAM_MAP = fits.ImageHDU(data=LOGLAM_MAP,header=header,name='LOGLAM_MAP')
                tbhdu = fits.BinTableHDU.from_columns(cols,header=header)

                hdulist = fits.HDUList([hdu_DELTA,hdu_iv,hdu_LOGLAM_MAP,tbhdu])
                hdulist.writeto(save_location + '/' + 'node_%s_nside_%d_pix_%d.fits' % (node, N_side, n))

            else:
                #Some kind of error and option to put in a new format code?
                print('Try another format!')

        else:
            print('No objects in pixel %d.' % n)

    return THING_ID, pixel_ID
