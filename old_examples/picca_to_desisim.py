import os
import numpy as np
import fitsio
import healpy as hp
import desitarget.mock.io as mockio
import lya_mock_functions as mock

# identify output file we want to plot
filename='../example_data/delta_picca/z1.85/z1.85_N1000_node_015_nside_4_pix_10.fits'
h = fitsio.FITS(filename)

# get information about quasars (TYPE,RA,DEC,Z_COSMO,DZ_RSD)
catalog = h[3].read()
# get arraw with redshift in each cell of grid
loglam = h[2].read()
# Get deltas (fluctuation around mean density)
delta = h[0].read()
z_qso = catalog['Z']
print('# initial quasars =',len(z_qso))
# keep only quasars with z>2.0
highz=z_qso>2.0
z_qso = z_qso[highz]
catalog = catalog[highz]
print('original shape',delta.shape)
delta = delta[:,highz]
print(np.min(z_qso),'< z_qso <',np.max(z_qso))
print('# high-z quasars =',len(z_qso))
# keep only quasars with good Dec
bad_dec = np.isnan(catalog['DEC']) | (catalog['DEC'] < -90.0) | (catalog['DEC'] > 90.0)
catalog = catalog[np.invert(bad_dec)]
z_qso = z_qso[np.invert(bad_dec)]
delta = delta[:,np.invert(bad_dec)]
Nq=len(z_qso)
print('# good quasars =',len(z_qso))
wave=np.power(10.0,loglam)
z = wave/1215.67-1.0
Nz=len(z)
print('full wavelength shape',delta.shape)
# we will only write pixels with wavelength in DESI spectrograph
in_desi=wave>3550.0
z = z[in_desi]
wave = wave[in_desi]
delta = delta[in_desi]
print('DESI shape',delta.shape)

# identify HEALPix pixels for our quasars
nside=8
Npix=12*nside*nside
pixels=hp.ang2pix(nside,(catalog['DEC']+90.0)/180.0*np.pi,catalog['RA']/180.0*np.pi)

transmission_base_dir='test_dir'

for pix in range(Npix):
    # get quasars in HEALPix pixel
    in_pix=(pixels==pix)
    # if none, move on
    if len(catalog[in_pix]) == 0:
        continue
    # select relevant quasars
    qso_in_pix = catalog[in_pix]
    delta_in_pix = delta[:,in_pix]
    z_in_pix = z_qso[in_pix]
    N_in_pix=len(qso_in_pix)
    print('useful pixel',pix,N_in_pix)
    # open file to write
    #filename='test_'+str(pix)+'.fits'
    dirname=mockio.get_healpix_dir(nside, pix, basedir=transmission_base_dir)
    #print('dirname',dirname)
    os.makedirs(dirname, exist_ok=True)
    filename=mockio.findfile('transmission', nside, pix, transmission_base_dir)
    #print('filename',filename)
    fits = fitsio.FITS(filename,'rw',clobber=True)

    for i,qso in enumerate(qso_in_pix):
        # Convert density to flux
        tau = mock.get_tau(z,1+delta_in_pix[:,i])
        flux = np.exp(-tau)
        # only add absorption in the forest
        no_forest = (z > z_in_pix[i])
        flux[no_forest]=1.0
        data = {}
        data['LAMBDA']=wave
        data['FLUX']=flux
        head = {}
        head['ZQSO']=z_in_pix[i]
        head['RA']=qso['RA']
        head['DEC']=qso['DEC']
        head['MAG_G']=22
        fits.write(data,header=head)
    fits.close()
