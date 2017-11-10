import numpy as np
from astropy.io import fits
from astropy.table import Table
import healpy as hp

#Program to convert a .fits file of HEALPix pixels into several .fits files of smaller HEALPix pixels.
#Takes files of the format outputted by CoLoRe and converts into the same format as outputted by ForetFusion.

#Function to split an individual .fits file into per-pixel .fits files.
#Takes the desried Nside, original .fits file filename and desired save location as inputs.
#Saves the per-pixel .fits files at the specified save location, and outputs THING_ID and pixel_ID lists.
def split_file(Nside,original_filename,save_location):
    #Open the original .fits file, and extract the coordinates of the objects.
    initial = fits.open(original_filename)
    RA = initial[1].data['RA']
    DEC = initial[1].data['DEC']

    #Convert the coordinates into new pixel identifier numbers, according to the Nside specified.
    pixel_ID=np.zeros([1,len(RA)])
    for i in range(len(RA)):
        #Can we just import healpy on cori? May need to install
        pixel_ID[0,i] = hp.pixelfunc.ang2pix(Nside,(np.pi/180.0)*(90.0-DEC[i]),(np.pi/180.0)*RA[i])

    #Set up a pixel_ID list to map between objects and their pixel.
    pixel_ID = pixel_ID.reshape(-1)

    #Set up a list of pixels represented in the original .fits file.
    pixel_list = list(set(pixel_ID))

    #NEEDS UPDATING to actual LOGLAM_MAP
    LOGLAM_MAP = np.zeros(initial[2].data[:].shape[1])
    #NEEDS UPDATING to actual THING_ID
    THING_ID = range(initial[2].data[:].shape[0])

    #Extract the density from the CoLoRe output file.
    DENSITY = np.transpose(initial[2].data[:])
    print(DENSITY.shape)

    #For each pixel represented in the original .fits file, make a new file.
    for n in pixel_list:

        #Set up the per-pixel THING_ID and DENSITY structures to go into the new per-pixel .fits file
        pixel_THING_ID = []
        pixel_DENSITY = np.zeros([DENSITY.shape[0],1])

        #For each object, assess if they are in the current pixel (n).
        for i in range(pixel_ID.shape[0]):

            #If so, take the data from the original and add it to the per-pixel output.
            if pixel_ID[i]==n:
                pixel_THING_ID = pixel_THING_ID+[THING_ID[i]]
                pixel_DENSITY = np.concatenate((pixel_DENSITY,np.transpose(np.asarray([DENSITY[:,i]]))),axis=1)


        #Remove the first column of the DENSITY array, which is blank.
        pixel_DENSITY = np.delete(pixel_DENSITY,0,1)


        if len(pixel_THING_ID)!=0:
            #Create hdus from the data arrays
            hdu_THING_ID = fits.ImageHDU(data=np.asarray(pixel_THING_ID),header=None,name='THING_ID')
            hdu_LOGLAM_MAP = fits.ImageHDU(data=LOGLAM_MAP,header=None,name='LOGLAM_MAP')
            hdu_DENSITY = fits.ImageHDU(data=pixel_DENSITY,header=None,name='DENSITY')

            #Create a primary HDU.
            prihdr = fits.Header()
            prihdr['COMMENT'] = ''
            prihdu = fits.PrimaryHDU(header=prihdr)

            hdulist = fits.HDUList([prihdu, hdu_THING_ID, hdu_LOGLAM_MAP, hdu_DENSITY])
            hdulist.writeto(save_location + 'pix_%d.fits' % n)

        else:
            print('No objects in pixel %d.' % n)

    return THING_ID, pixel_ID


#Function to take multiple input files, split them each with 'split_file' and create a master.fits file.
def split_pixels(Nside,files,save_location):

    #Determine number of files to be split.
    N_files = len(files)

    #Set up the structures for the master file.
    master_pixel_ID=[]
    master_THING_ID=[]

    #Run through each file and split it into pixel fits files, adding each pixel's THING_IDs and pixel_IDs to the master.
    for i in range(N_files):
        THING_ID, pixel_ID = split_file(Nside,files[i],save_location)
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
