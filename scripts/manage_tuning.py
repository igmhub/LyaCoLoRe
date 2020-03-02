import numpy as np
from astropy import fits

z_values = np.linspace(2.0,3.375,12)

#inputs: n, k1, number of processes, overwrite or not?,filename if not,add or not?, number of pixels, filename of old ones

results_list = []

#Load previous data
h = fits.open('input_files/tuning.fits')
tuning_data = h[1].data
tuning_data = np.sort(tuning_data,order=['z'])

if not n:
    n = h[1].header['n']
if not k1:
    k1 = h[1].header['k1']

for z_value in z_values:
    a_start = np.interp(z_value,tuning_data['z'],tuning_data['alpha'])
    b_start = tuning_data['beta'][i]
    sG_start = tuning_data['sigma_G'][i]
    print(z_value,alpha_start,beta_start,sigma_G_start)
    #minuit = iminuit_tuning.tune(z_value,alpha_start,beta_start,sigma_G_start,N_processes,N_pixels)
    #results_list += [(z_value,minuit.values['alpha'],minuit.values['beta'],minuit.values['sigma_G'])]

dtype = [('z', 'f4'), ('alpha', 'f4'), ('beta', 'f4'), ('sigma_G', 'f4')]
results = np.array(results_list,dtype=dtype)

if combine_results:
    to_save = np.concatenate((tuning_data,results))
    to_save = np.sort(to_save,order=['z'])
else:
    to_save = np.sort(results,order=['z'])

header = fits.Header()
header['n'] = n
header['k1'] = k1

prihdr = fits.Header()
prihdu = fits.PrimaryHDU(header=prihdr)
cols_DATA = fits.ColDefs(to_save)
hdu_DATA = fits.BinTableHDU.from_columns(cols_DATA,header=header,name='DATA')

hdulist = fits.HDUList([prihdu, hdu_DATA])
hdulist.writeto('input_files/tuning.fits',overwrite=True))
hdulist.close()
