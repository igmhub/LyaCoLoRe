import numpy as np
from astropy.io import fits

basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_revamp/process_output_4096_32/'
basedir = '/Users/jfarr/Projects/repixelise/test_output/test_multi/'
quantities = ['gaussian','density','flux']
N_side = 8
pixels = [3,10,11,21,22,23,36,37,38,39,56,57,58,80,81,108,428,459,460,491,492,493,522,523,524,525,555,556,557,587,588,620]

for pixel in pixels:
    pixel_100 = pixel//100
    N_cells = []
    for quantity in quantities:
        h = fits.open(basedir+'{}/{}/picca-{}-{}-{}.fits'.format(pixel_100,pixel,quantity,N_side,pixel))
        final_cell_row = h[0].data[-1,:]
        final_IVAR_row = h[1].data[-1,:]
        N_cells += [h[0].shape[0]]
        h.close()
    if sum(N_cells) != 571*3:
        print('\n!!!!!')
        print(pixel,N_cells,end='\r')
        print('\n!!!!!')
    else:
        print(pixel,N_cells,end='\r')

    if sum(final_IVAR_row) == 0:
        print('final cells of {} have zero weight'.format(pixel))


print(' ')
