import os
import numpy as np
import time
from astropy.io import fits
import healpy as hp
import desitarget.mock.io as mockio
import lya_mock_functions as mock

# identify output file we want to plot
filename='/global/cscratch1/sd/jfarr/LyaSkewers/healpix_4096_32/per_pixel/nside_4_pix_1.fits'
#filename='/global/cscratch1/sd/jfarr/LyaSkewers/healpix_4096_32/per_pixel/nside_4_pix_0.fits'
h = fits.open(filename)

# get information about quasars (TYPE,RA,DEC,Z_COSMO,DZ_RSD)
catalog = h[3].data
# get arraw with redshift in each cell of grid
loglam = h[2].data
# get deltas (fluctuation around mean density) 
delta = h[0].data
print('original delta shape',delta.shape)

# we will only write pixels with wavelength in DESI spectrograph
wave=np.power(10.0,loglam)
z = wave/1215.67-1.0
in_desi = wave>3550.0
in_desi = np.where(in_desi)[0]
wave = wave[in_desi]
z = z[in_desi]
t0=time.time()
delta = delta[in_desi]
t1=time.time()
print('desi delta shape',delta.shape,'time =',t1-t0)

# keep only quasars with z>2.0
z_qso = catalog['Z']
print('# initial quasars =',len(z_qso))
print(np.min(z_qso),'< z_qso <',np.max(z_qso))
highz = z_qso>2.0
highz = np.where(highz)[0]
z_qso = z_qso[highz]
print('# high-z quasars =',len(z_qso))
catalog = catalog[highz]
t0=time.time()
delta = delta[:,highz]
t1=time.time()
print('highz delta shape',delta.shape,'time =',t1-t0)

# keep only quasars with good Dec
bad_dec = np.isnan(catalog['DEC']) | (catalog['DEC'] < -90.0) | (catalog['DEC'] > 90.0)
good_dec = np.invert(bad_dec)
good_dec = np.where(good_dec)[0]
catalog = catalog[good_dec]
z_qso = z_qso[good_dec]
Nq=len(z_qso)
print('# good quasars =',len(z_qso))
t0=time.time()
delta = delta[:,good_dec]
t1=time.time()
print('final delta shape',delta.shape,'time =',t1-t0)

# identify HEALPix pixels for our quasars
nside=16
Npix=12*nside*nside
pixels=hp.ang2pix(nside,(catalog['DEC']+90.0)/180.0*np.pi,catalog['RA']/180.0*np.pi)

transmission_base_dir='test_dir_'+str(nside)
for pix in range(Npix):
    # get quasars in HEALPix pixel
    in_pix=(pixels==pix)
    # if none, move on
    if len(catalog[in_pix]) == 0:
        continue
    # select relevant quasars
    qso_in_pix = catalog[in_pix]
    delta_in_pix = delta[:,in_pix]
    N_in_pix=len(qso_in_pix)
    print('useful pixel',pix,N_in_pix)

    #Add a couple of headers to the file.
    header = fits.Header()
    header['NSIDE'] = nside
    header['NQSO'] = N_in_pix
    header['PIX'] = int(pix)
    prim_hdu = fits.PrimaryHDU(header=header)

    #meta-data
    ra=qso_in_pix['RA']
    dec=qso_in_pix['DEC']
    zq=qso_in_pix['Z']
    mockid=qso_in_pix['THING_ID']
    
    #Construct a table for the meta-data hdu
    col_ra = fits.Column(name='RA', array=ra, format='E')
    col_dec = fits.Column(name='DEC', array=dec, format='E')
    col_zq = fits.Column(name='Z', array=zq, format='E')
    col_id = fits.Column(name='MOCKID', array=mockid, format='A10')
    meta_hdu = fits.BinTableHDU.from_columns([col_ra, col_dec, col_zq, col_id],name='METADATA')
    #meta_hdu.writeto('test_'+str(pix)+'.fits')
        
    flux = np.empty_like(delta_in_pix)
    for i in range(N_in_pix):
        # Convert density to flux
        tau = mock.get_tau(z,1+delta_in_pix[:,i])
        toflux = np.exp(-tau)
        # only add absorption in the forest 
        no_forest = (z > z_qso[i])
        toflux[no_forest]=1.0
        flux[:,i]=toflux
    
    wave_hdu = fits.ImageHDU(data=wave,name='WAVELENGTH')
    flux_hdu = fits.ImageHDU(data=flux,name='TRANSMISSION')
    hdulist = fits.HDUList([prim_hdu,meta_hdu,wave_hdu,flux_hdu])

    # open file to write
    #filename='test_'+str(pix)+'.fits'
    dirname=mockio.get_healpix_dir(nside, pix, basedir=transmission_base_dir)
    #print('dirname',dirname)
    os.makedirs(dirname, exist_ok=True)
    filename=mockio.findfile('transmission', nside, pix, transmission_base_dir)
    #print('filename',filename)
    #fits = fitsio.FITS(filename,'rw',clobber=True)
    hdulist.writeto(filename)

master_file=transmission_base_dir+"/master.fits"
print('master file',master_file)

#Add a couple of headers to the file.
header = fits.Header()
header['NSIDE'] = nside
header['NQSO'] = Nq

#Construct a table for the meta-data hdu
col_ra = fits.Column(name='RA', array=catalog['RA'], format='E')
col_dec = fits.Column(name='DEC', array=catalog['DEC'], format='E')
col_zq = fits.Column(name='Z', array=catalog['Z'], format='E')
col_id = fits.Column(name='MOCKID', array=catalog['THING_ID'], format='A10')
col_pix = fits.Column(name='PIXNUM', array=pixels, format='J')
meta_hdu = fits.BinTableHDU.from_columns([col_ra, col_dec, col_zq, col_id, col_pix],name='METADATA',header=header)
meta_hdu.writeto(master_file)
