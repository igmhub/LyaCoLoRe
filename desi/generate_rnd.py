import numpy as np
import astropy.table
from scipy.interpolate import interp1d
import healpy as hp
from astropy.io import fits
import warnings

from lyacolore import catalog,utils

#Set up options
factor = 10.0
out_path = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v9/v9.0.0_full/master_randoms.fits'
#out_path = '/global/projecta/projectdirs/desi/mocks/lya_forest/develop/london/v8.0/v8.0.0/master_randoms.fits'
#out_path = '/project/projectdirs/desi/users/jfarr/LyaCoLoRe_paper/data/LyaCoLoRe_output/v9.0.9/master_randoms.fits'
method = 'from_catalog'
catalog_path = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v9/v9.0.0_full/master.fits'
#catalog_path = '/global/projecta/projectdirs/desi/mocks/lya_forest/develop/london/v8.0/v8.0.0/master.fits'
#catalog_path = '/project/projectdirs/desi/users/jfarr/LyaCoLoRe_paper/data/LyaCoLoRe_output/v9.0.9/master.fits'
footprint = 'full_sky'
nz_filename = 'input_files/Nz_qso_130618_2_colore1_hZs.txt'
min_cat_z = 1.8
max_cat_z = 3.79
overwrite = True
N_side = 16
start_MOCKID_rnd = 10**10

def generate_rnd(factor=3, out_path= None, method='use_catalog', catalog_path=None, footprint=None, nz_filename='input_files/Nz_qso_130618_2_colore1_hZs.txt', min_cat_z=1.8, max_cat_z=4.0, overwrite=False, N_side=16,start_MOCKID_rnd=10**10):
    """
    Routine to generate a random catalog in 3D following
    certain N(z) distribution

    Args:
    ----
    factor: Size of the generated catalog (before masking)
    out_path: Name of output file where randoms will be saved (default: None)
    method: Method to generate the random catalog (default: 'random_choice')
    """

    #Set up vectors from n(z) file.
    N_vec = 500
    nz_file = astropy.table.Table.read(nz_filename,format='ascii')
    zedges = np.linspace(min_cat_z,max_cat_z,N_vec)
    zvec = (zedges[1:] + zedges[:-1])/2.
    dz = zvec[1] - zvec[0]

    #Get the total number of QSOs.
    ntot = 4*180**2/np.pi*np.sum(nz_file['col2'])*(nz_file['col1'][1]-nz_file['col1'][0])
    ntot = int(factor*ntot)

    #Generate random redshifts by one of 3 different methods:
    if method=='from_catalog':
        #Method 1: choose from the QSO catalog and add small deviations.
        #Make a z vector with edges extended from min and max cat values.
        #This is to take into account RSDs in QSOs near the edge of the range.
        pc_extra = 0.05
        extra = (max_cat_z - min_cat_z) * pc_extra
        zedges = np.linspace(min_cat_z-extra,max_cat_z+extra,N_vec*(1+2*extra))
        zvec = (zedges[1:] + zedges[:-1])/2.
        dz = zvec[1] - zvec[0]

        #Get catalog data and calculate dndz from it.
        if catalog_path is None:
            raise ValueError('Needed a path to read the catalog')
        tab = fits.open(catalog_path)['CATALOG'].data
        z_master_RSD = tab['Z_QSO_RSD']
        dndz_RSD,_ = np.histogram(z_master_RSD,bins=zedges)

        #Get ntot
        ntot = factor*tab.shape[0]

        #Turn this into a cdf and draw redshifts from it
        cdf_RSD = np.cumsum(dndz_RSD)/np.sum(dndz_RSD)
        cdf_RSD_i = np.concatenate([[0],cdf_RSD])
        icdf_RSD = interp1d(cdf_RSD_i,zedges,fill_value=(0.,1.),bounds_error=False)
        z_rnd = icdf_RSD(np.random.random(size=ntot)) + np.random.normal(size=ntot,scale=10**-6)

    elif method=='rnd_choice':
        #Method 2: choose from within the vectorisation of the dndz file.
        #This just seems like a less good version of 3?
        z_rnd = np.random.choice(zvec,p=dndz/np.sum(dndz),size=ntot) + 2*np.random.normal(size=ntot,scale=dz)

    elif method=='cdf':
        #Method 3: Calculate the cumulative distribution and interpolate.
        #Get dndz by interpolating input file..
        spl_z = interp1d(nz_file['col1'],nz_file['col2'],fill_value=0.)
        dndz = spl_z(zvec)

        #Generate redshifts without RSDs.
        cdf = np.cumsum(dndz)/np.sum(dndz)
        cdf_i = np.concatenate([[0],cdf])
        icdf = interp1d(cdf_i,zedges,fill_value=(0.,1.),bounds_error=False)
        z_rnd = icdf(np.random.random(size=ntot))

        #Measure sigma_RSD from the catalog and apply it as a Gaussian shift.
        if catalog_path is None:
            raise ValueError('Needed a path to read the catalog')
        tab = fits.open(catalog_path)['CATALOG'].data
        z_master_NO_RSD = tab['Z_QSO_NO_RSD']
        z_master_RSD = tab['Z_QSO_RSD']
        dz_rsd_master = z_master_RSD - z_master_NO_RSD
        sigma_rsd_master = np.std(dz_rsd_master)
        z_rnd += np.random.normal(size=ntot,scale=sigma_rsd_master)

    #Assign random positions on the sky to the QSOs.
    ra_rnd = 360.*np.random.random(size=len(z_rnd))
    cth_rnd = -1+2.*np.random.random(size=len(z_rnd))
    dec_rnd = np.arcsin(cth_rnd)*180/np.pi

    #Assign MOCKIDs to the QSOs, checking that there's no overlap with those in
    #the master file.
    if catalog_path is None:
        raise ValueError('Needed a path to read the catalog')
    tab = fits.open(catalog_path)['CATALOG'].data
    max_cat_MOCKID = np.max(tab['MOCKID'])
    while max_cat_MOCKID > start_MOCKID_rnd:
        warnings.warn('Start value of randoms\' MOCKIDs is not high enough: increasing from {} to {}'.format(start_MOCKID_rnd,10*start_MOCKID_rnd))
        start_MOCKID_rnd *= 10
    MOCKID_rnd = (np.array(list(range(ntot))) + start_MOCKID_rnd).astype(int)
    if np.max(MOCKID_rnd) > (2**63 - 1):
        raise ValueError('Max MOCKID exceeds max integer allowed by FITS.')

    #Filter the QSOs according to the input footprint.
    QSO_filter = utils.make_QSO_filter(footprint)
    good = QSO_filter(ra_rnd,dec_rnd)
    ra_rnd = ra_rnd[good]
    dec_rnd = dec_rnd[good]
    MOCKID_rnd = MOCKID_rnd[good]
    z_rnd = z_rnd[good]
    pix_rnd = utils.make_pixel_ID(N_side,ra_rnd,dec_rnd)

    #Write the catalog to file.
    dtype = [('RA', 'd'), ('DEC', 'd'), ('Z', 'd'), ('PIXNUM', int), ('MOCKID', int)]
    ID_data = np.array(list(zip(ra_rnd,dec_rnd,z_rnd,pix_rnd,MOCKID_rnd)),dtype=dtype)
    if out_path is not None:
        catalog.write_ID(out_path,N_side,ID_data,overwrite=overwrite)

    return

# Execute
generate_rnd(factor=factor,out_path=out_path,method=method,catalog_path=catalog_path,footprint=footprint,nz_filename=nz_filename,min_cat_z=min_cat_z,max_cat_z=max_cat_z,overwrite=overwrite,N_side=N_side,start_MOCKID_rnd=start_MOCKID_rnd)
