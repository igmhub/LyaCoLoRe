import numpy as np
import astropy.table
from scipy.interpolate import interp1d
import healpy as hp
from astropy.io import fits

from lyacolore import catalog,utils

#Set up options
factor = 10
out_path = '/global/projecta/projectdirs/desi/mocks/lya_forest/develop/london/v7.3/v7.3.0/master_randoms.fits'
method = 'cdf'
catalog_path = '/global/projecta/projectdirs/desi/mocks/lya_forest/develop/london/v7.3/v7.3.0/master.fits'
footprint = 'desi_pixel_plus'
nz_filename = 'input_files/Nz_qso_130618_2_colore1_hZs.txt'
min_cat_z = 1.8
max_cat_z = 3.79
overwrite = True
N_side = 16

def generate_rnd(factor=3, out_path= None, method='use_catalog', catalog_path=None, footprint=None, nz_filename='input_files/Nz_qso_130618_2_colore1_hZs.txt', min_cat_z=1.8, max_cat_z=3.79, overwrite=False, N_side=16):
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

    #Get the input n(z) file and set up vectors.
    nz_file = astropy.table.Table.read(nz_filename,format='ascii')
    zvec = np.linspace(np.min(nz_file['col1']),np.max(nz_file['col1']),5000)
    spl_z = interp1d(nz_file['col1'],nz_file['col2'])
    dndz = spl_z(zvec)

    #Fiter by the minimum z.
    ind = (zvec>=min_cat_z) * (zvec<=max_cat_z)
    zvec = zvec[ind]
    dndz = dndz[ind]
    dz = zvec[1] - zvec[0]

    #Get the total number of QSOs.
    ntot = 4*180**2/np.pi*np.sum(nz_file['col2'])*(nz_file['col1'][1]-nz_file['col1'][0])
    ntot = int(factor*ntot)

    #Generate random redshifts by one of 3 different methods:
    if method=='use_catalog':
        #Method 1: choose from the QSO catalog and add small deviations.
        if catalog_path is None:
            raise ValueError('Needed a path to read the catalog')
        tab = fits.open(catalog_path)[1].data
        z_rnd = np.random.choice(tab['Z_QSO_NO_RSD'], size=ntot)+1e-6*np.random.normal(size=ntot)

    elif method=='rnd_choice':
        #Method 2: choose from within the vectorisation of the dndz file.
        z_rnd = np.random.choice(zvec,p=dndz/np.sum(dndz),size=ntot) + 2*np.random.normal(size=ntot,scale=dz)

    elif method=='cdf':
        #Method 3: Calculate the cumulative distribution and interpolate.
        cdf = np.cumsum(dndz)/np.sum(dndz)
        icdf = interp1d(cdf,zvec,fill_value='extrapolate',bounds_error=False)
        z_rnd = icdf(np.random.random(size=ntot))

    #Assign random positions on the sky to the QSOs.
    ra_rnd = 360.*np.random.random(size=len(z_rnd))
    cth_rnd = -1+2.*np.random.random(size=len(z_rnd))
    dec_rnd = np.arcsin(cth_rnd)*180/np.pi

    #Filter the QSOs according to the input footprint.
    QSO_filter = utils.make_QSO_filter(footprint)
    good = QSO_filter(ra_rnd,dec_rnd)
    ra_rnd = ra_rnd[good]
    dec_rnd = dec_rnd[good]
    z_rnd = z_rnd[good]
    pix_rnd = utils.make_pixel_ID(N_side,ra_rnd,dec_rnd)

    #Write the catalog to file.
    dtype = [('RA', 'd'), ('DEC', 'd'), ('Z', 'd'), ('PIXNUM', int)]
    ID_data = np.array(list(zip(ra_rnd,dec_rnd,z_rnd,pix_rnd)),dtype=dtype)
    if out_path is not None:
        catalog.write_ID(out_path,N_side,ID_data,overwrite=overwrite)

    return

# Execute
generate_rnd(factor=factor,out_path=out_path,method=method,catalog_path=catalog_path,footprint=footprint,nz_filename=nz_filename,min_cat_z=min_cat_z,max_cat_z=max_cat_z,overwrite=overwrite,N_side=N_side)
