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

    #Set up the structures for the master file.
    master_pixel_ID=[]
    master_THING_ID=[]

    #Run through each file and split it into pixel fits files, adding each pixel's THING_IDs and pixel_IDs to the master.
    for i in range(N_files):
        THING_ID, pixel_ID = split_file(N_side,files[i],file_numbers[i],save_location,output_format)
        THING_ID = list(THING_ID)
        pixel_ID = list(pixel_ID)
        master_THING_ID += (THING_ID)
        master_pixel_ID += (pixel_ID)

    #Sort the THING_IDs and pixel_IDs into the right order, as specified in ForetFusion.
    master_ID = np.array([master_THING_ID,master_pixel_ID])
    master_ID_T = np.transpose(master_ID)
    master_ID_T_sort = np.sort(master_ID_T,axis=0)
    master_ID_T_sort = np.sort(master_ID_T,axis=-1)
    master_ID_sort = np.transpose(master_ID_T_sort)

    #Take the output from each splitting and make it into a master file
    master = Table([master_ID_sort[0,:],master_ID_sort[1,:]], names=('THING_ID','PIX'))
    master.write(save_location + 'master.fits', format='fits')

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
    DENSITY = initial[2].data[:]

    N_cells = z.shape[0]
    N_qso = z_qso.shape[0]
    N_pix = 12*N_side**2

    iv = np.zeros((N_qso,N_cells))
    PLATE = np.zeros(N_qso)
    MJD = np.zeros(N_qso)
    FIBER = np.zeros(N_qso)



    #Convert the coordinates into new pixel identifier numbers, according to the N_side specified.
    pixel_ID=np.zeros([1,len(RA)])
    for i in range(len(RA)):
        #Can we just import healpy on cori? May need to install
        pixel_ID[0,i] = hp.pixelfunc.ang2pix(N_side,(np.pi/180.0)*(90.0-DEC[i]),(np.pi/180.0)*RA[i])

    #Set up a pixel_ID list to map between objects and their pixel.
    pixel_ID = pixel_ID.reshape(-1)

    #Set up a list of pixels represented in the original .fits file.
    pixel_list = list(set(pixel_ID))

    #NEEDS UPDATING to actual LOGLAM_MAP
    LOGLAM_MAP = 1215.67*(1+z)

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



    #For each pixel represented in the original .fits file, make a new file.
    for n in pixel_list:

        #Progress check aide.
        print('Working on pixel %d of %d...' % (n+1,N_pix))

        #Set up the per-pixel structures to go into the new per-pixel .fits file
        pixel_DENSITY = np.zeros([1,N_cells])
        pixel_iv = np.zeros([1,N_cells])
        pixel_RA = []
        pixel_DEC = []
        pixel_z_qso = []
        pixel_PLATE = []
        pixel_MJD = []
        pixel_FIBER = []
        pixel_THING_ID = []

        for i in range(pixel_ID.shape[0]):
            #For each object, assess if they are in the current pixel (n).
            #If so, take the data from the original and add it to the per-pixel output.
            if pixel_ID[i]==n:
                pixel_DENSITY = np.concatenate((pixel_DENSITY,np.asarray([DENSITY[i,:]])),axis=0)
                pixel_iv = np.concatenate((pixel_iv,np.asarray([iv[i,:]])),axis=0)
                pixel_RA += [RA[i]]
                pixel_DEC += [DEC[i]]
                pixel_z_qso += [z_qso[i]]
                pixel_PLATE += [PLATE[i]]
                pixel_MJD += [MJD[i]]
                pixel_FIBER += [FIBER[i]]
                pixel_THING_ID += [THING_ID[i]]

        #Remove the first row of the DENSITY array, which is blank.
        pixel_DENSITY = np.delete(pixel_DENSITY,0,0)
        pixel_iv = np.delete(pixel_iv,0,0)


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

                #Create hdus from the data arrays
                hdu_DENSITY = fits.PrimaryHDU(data=pixel_DENSITY)
                hdu_iv = fits.ImageHDU(data=pixel_iv,header=None,name='IV')
                hdu_LOGLAM_MAP = fits.ImageHDU(data=LOGLAM_MAP,header=None,name='LOGLAM_MAP')
                tbhdu = fits.BinTableHDU.from_columns(cols)

                hdulist = fits.HDUList([hdu_DENSITY,hdu_iv,hdu_LOGLAM_MAP,tbhdu])
                hdulist.writeto(save_location + 'pix_%d.fits' % n)
            else:
                #Some kind of error and option to put in a new format code?
                print('Try another format!')
        else:
            print('No objects in pixel %d.' % n)

    return THING_ID, pixel_ID
