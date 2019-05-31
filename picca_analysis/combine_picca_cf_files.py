#!/usr/bin/env python
  
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from multiprocessing import Pool
import multiprocessing
import argparse
import glob

################################################################################

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--in-dir', type = str, default = None, required=True,
                    help = 'input data directory')

parser.add_argument('--out', type = str, default = None, required=True,
                    help = 'output filepath')

parser.add_argument('--files', type = str, default = None, required=False,
                    help = 'list of files to combine', nargs='*')

parser.add_argument('--overwrite', action="store_true", default = False, required=False,
                    help = 'overwrite existing files')

################################################################################

print('setup arguments from parser')

args = parser.parse_args()

base_in_dir = args.in_dir
out_filepath = args.out
if args.files:
    file_list = [base_in_dir+f for f in args.files]
else:
    file_list = glob.glob(base_in_dir+'cf*.fits.gz')
overwrite = args.overwrite

################################################################################

#Open all of the files and add the data to a common structure.
h = fits.open(file_list[0])
prihdu_header = h[0].header
hdu_1_data = h[1].data
hdu_1_header = h[1].header
hdu_2_data = h[2].data
hdu_2_header = h[2].header
h.close()

for f in file_list[1:]:
    h = fits.open(f)

    #Maybe ought to check that the RP and RT bins are the same here?
    hdu_1_data['NB'] += h[1].data['NB']
    hdu_2_data = np.concatenate((hdu_2_data,h[2].data))
    hdu_2_header['NAXIS2'] += 1

    h.close()

#Make the data into suitable HDUs.
PRIHDU = fits.PrimaryHDU(header=prihdu_header)
HDU_1 = fits.BinTableHDU(data=hdu_1_data,header=hdu_1_header,name='ATTRI')
HDU_2 = fits.BinTableHDU(data=hdu_2_data,header=hdu_2_header,name='COR')

#Combine the HDUs into and HDUlist and save as a new file. Close the HDUlist.
hdulist = fits.HDUList([PRIHDU, HDU_1, HDU_2])
hdulist.writeto(out_filepath,overwrite=overwrite)
hdulist.close()

