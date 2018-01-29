import numpy as np
from astropy.io import fits
import sys

"""
Function to take a master file outputted by 'process_colore_multi', and convert
it into a similar file which can be read by picca's "do_xcf" cross-correlation
function.
"""

original_filename = sys.argv[1]

#Can use 'Z_QSO' or 'Z_QSO_NO_RSD'
RSD_option = sys.argv[2]

#lambda_min = 3550.0
#lya = 1215.67
#z_min = lambda_min/lya - 1.0

#Set things up to include RSDs or not.
if RSD_option == 'RSD':
    Z_option = 'Z_QSO'
elif RSD_option == 'NO_RSD':
    Z_option = 'Z_QSO_NO_RSD'
else:
    error('RSD_option (second argument) not recognised. Please enter \'RSD\' or \'NO_RSD\'. ')

new_filename = original_filename[:-5] + '_picca_{}.fits'.format(RSD_option)

#Import data from the original master file.
#This relies on the master file being of the format produced by 'process_colore_multi'.
original = fits.open(original_filename)
header = original[1].header
RA = original[1].data['RA']
DEC = original[1].data['DEC']
THING_ID = original[1].data['MOCKID']
Z = original[1].data[Z_option]

N_qso = RA.shape[0]

#Produce the extra 'data' that we need for running picca do_xcf.py.
#Id does not seem to actually use this data, but will crash if columns containing it do not exist.
MJD = np.zeros(N_qso)
FID = np.zeros(N_qso)
PLATE = THING_ID

dtype = [('RA','>f4'),('DEC','>f4'),('Z','>f4'),('THING_ID',int),('MJD','>f4'),('FIBERID',int),('PLATE',int)]
new_data = np.array(list(zip(RA,DEC,Z,THING_ID,MJD,FID,PLATE)),dtype=dtype)

#Create a new master file, with the same filename concatenated with '_picca_' and the RSD option chosen.
prihdr = fits.Header()
prihdu = fits.PrimaryHDU(header=prihdr)

bt = fits.BinTableHDU.from_columns(new_data,header=header)

hdulist = fits.HDUList([prihdu,bt])
hdulist.writeto(new_filename)
hdulist.close()
