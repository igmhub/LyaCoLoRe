import numpy as np
from desimodel.footprint import is_point_in_desi
from desimodel.io import load_tiles
import astropy.table
from scipy.interpolate import interp1d
import healpy as hp
import fitsio
from desimodel import footprint
from desimodel.io import load_pixweight
import os
def generate_rnd(factor=3, out_path= None, method='use_catalog', catalog_path= None):
    """
    Routine to generate a random catalog in 3D following
    certain N(z) distribution
    
    Args:
    ----
    rad: z-values of the data catalog
    factor: Size of the generated catalog (before masking)
    out_path: Name of output file where randoms will be saved (default: None)
    method: Method to generate the random catalog (default: 'random_choice')
    """
    #Creating random that follows N(z)
    footprint_healpix_nside=256
    tiles = load_tiles()
    nz_file = astropy.table.Table.read('example_data/Nz_qso_2_highZ.txt',format='ascii')
    zvec = np.linspace(np.min(nz_file['col1']),np.max(nz_file['col1']),5000)
    spl_z = interp1d(nz_file['col1'],nz_file['col2'])
    dndz = spl_z(zvec)
    ntot = 4*180**2/np.pi*np.sum(nz_file['col2'])*(nz_file['col1'][1]-nz_file['col1'][0])
    ntot = int(factor*ntot)
    if method=='use_catalog':
        if catalog_path is None:
            raise ValueError('Needed a path to read the catalog')
        tab = fitsio.read(catalog_path)
        z_rnd = np.random.choice(tab['Z_QSO_RSD'], size=ntot)+1e-6*np.random.normal(size=ntot)  
    if method=='rnd_choice':
        z_rnd = np.random.choice(zvec[zvec>1.8],p=dndz[zvec>1.8]/np.sum(dndz[zvec>1.8]),size=ntot)+2*(z[1]-z[0])*np.random.normal(size=ntot)
    if method=='cdf':
        cdf = np.cumsum(dndz[zvec>=1.8])/np.sum(dndz[zvec>=1.8])
        icdf = interp1d(cdf,zvec[zvec>=1.8],fill_value='extrapolate',bounds_error=False)
        z_rnd = icdf(np.random.random(size=ntot))
    #tab = fitsio.read('/global/projecta/projectdirs/desi/mocks/lya_forest/london/v4.0/quick-0.0/zcat_desi_drq.fits')
    ra_rnd = 360.*np.random.random(size=len(z_rnd))
    cth_rnd = -1+2.*np.random.random(size=len(z_rnd))
    dec_rnd = np.arcsin(cth_rnd)*180/np.pi
    #pixnums = hp.ang2pix(256,np.pi/2-tab['DEC']*np.pi/180, np.pi/180*tab['RA'])
    #pixnums2 = hp.ang2pix(256,np.pi/2-dec_rnd*np.pi/180, np.pi/180*ra_rnd) 
    footprint_filename=os.path.join(os.environ['DESIMODEL'],'data','footprint','desi-healpix-weights.fits')   
    pixmap=fitsio.read(footprint_filename)
    footprint_healpix_weight = load_pixweight(footprint_healpix_nside, pixmap=pixmap)
    footprint_healpix = footprint.radec2pix(footprint_healpix_nside, ra_rnd, dec_rnd)
    good = np.where(footprint_healpix_weight[footprint_healpix]>0.99)[0]
    #good = np.in1d(pixnums2,pixnums)
    #good = is_point_in_desi(tiles,ra_rnd,dec_rnd)
    ra_rnd = ra_rnd[good]
    dec_rnd = dec_rnd[good]
    z_rnd = z_rnd[good]
    if out_path is not None:
        tab_out = astropy.table.Table([ra_rnd,dec_rnd,z_rnd],names=('RA','DEC','Z'))
        tab_out.write(out_path,overwrite=True)
    return None
# Execute
generate_rnd(factor=10,out_path='/global/projecta/projectdirs/desi/mocks/lya_forest/london/v5.0.0/master_randoms.fits',method='use_catalog', catalog_path='/global/projecta/projectdirs/desi/mocks/lya_forest/london/v5.0.0/master.fits')
