import numpy as np
from scipy.interpolate import interp1d
import astropy.io.fits as fits

from lyacolore import DLA, utils

lya = utils.lya_rest

#Set up options
factor = 0.1
out_path = '/global/projecta/projectdirs/desi/mocks/lya_forest/develop/london/v7.3/v7.3.0/master_DLA_randoms.fits'
method = 'cdf'
DLA_catalog_path = '/global/projecta/projectdirs/desi/mocks/lya_forest/develop/london/v7.3/v7.3.0/master_DLA.fits'
QSO_catalog_path = '/global/projecta/projectdirs/desi/mocks/lya_forest/develop/london/v7.3/v7.3.0/master.fits'
footprint = 'desi_pixel_plus'
lambda_min = 3470.
lambda_max = 6550.
NHI_min = 17.2
NHI_max = 22.5
overwrite = True
N_side = 16
add_NHI = True

def generate_rnd(factor=3, out_path=None , DLA_catalog_path=None, QSO_catalog_path=None, footprint='desi_pixel_plus', lambda_min=3470., lambda_max=6550., NHI_min=17.2, NHI_max=22.5, overwrite=False, N_side=16, add_NHI=True, method='cdf'):
    """
    Routine to generate a random catalog in 3D following
    certain N(z) distribution

    Args:
    ----
    factor: Size of the generated catalog (before masking)
    out_path: Output path
    """

    #Generate a z vector and the dn/dz function.
    zmin = lambda_min/lya - 1
    zmax = lambda_max/lya - 1
    n_vec = 100
    zvec = np.linspace(zmin,zmax,n_vec)
    zedges = np.concatenate([[zvec[0]]-(zvec[1]-zvec[0])/2.,(zvec[1:]+zvec[:-1])*0.5,[zvec[-1]+(-zvec[-2]+zvec[-1])*0.5]]).ravel()
    dz = zvec[1] - zvec[0]
    dndz = DLA.dndz(zvec,NHI_min=NHI_min,NHI_max=NHI_max)

    #Use dndz to get the mean number of DLAs in each vector cell. Scale it by
    mean_n = dz * dndz
    mean_n *= factor

    #Get data about the QSO sample.
    h = fits.open(QSO_catalog_path)
    RA = h['CATALOG'].data['RA']
    DEC = h['CATALOG'].data['DEC']
    z_qso = h['CATALOG'].data['Z_QSO_NO_RSD']
    z_qso_rsd = h['CATALOG'].data['Z_QSO_RSD']
    pixnum = h['CATALOG'].data['PIXNUM']
    MOCKID = h['CATALOG'].data['MOCKID']
    n_qso = z_qso.shape[0]
    h.close()

    ntot = (np.sum(mean_n) * n_qso).astype('int')
    print('generating {} random DLAs...'.format(ntot))

    #Generate random redshifts for the DLAs.
    if method=='from_catalog':
        #Method 1: choose from the DLA catalog and add small deviations.
        #Make a z vector with edges extended from min and max cat values.
        #This is to take into account RSDs in QSOs near the edge of the range.
        pc_extra = 0.05
        extra = (max_cat_z - min_cat_z) * pc_extra
        zedges = np.linspace(min_cat_z-extra,max_cat_z+extra,N_vec*(1+2*extra))
        zvec = (zedges[1:] + zedges[:-1])/2.
        dz = zvec[1] - zvec[0]

        #Get catalog data and calculate dndz from it.
        if DLA_catalog_path is None:
            raise ValueError('Needed a path to read the catalog')
        tab = fits.open(DLA_catalog_path)['CATALOG'].data
        z_master_RSD = tab['Z_DLA_RSD']
        dndz_RSD,_ = np.histogram(z_master_RSD,bins=zedges)

        #Turn this into a cdf and draw redshifts from it
        cdf_RSD = np.cumsum(dndz_RSD)/np.sum(dndz_RSD)
        cdf_RSD_i = np.concatenate([[0],cdf_RSD])
        icdf_RSD = interp1d(cdf_RSD_i,zedges,fill_value=(0.,1.),bounds_error=False)
        z_rnd = icdf_RSD(np.random.random(size=ntot)) + np.random.normal(size=ntot,scale=10**-6)

    elif method=='cdf':
        #Method 2: Calculate the cumulative distribution and interpolate.
        #Get the input n(z) file and set up vectors.
        dndz = DLA.dndz(zvec,NHI_min=NHI_min,NHI_max=NHI_max)

        #Generate redshifts without RSDs.
        cdf = np.cumsum(dndz)/np.sum(dndz)
        cdf_i = np.concatenate([[0],cdf])
        icdf = interp1d(cdf_i,zedges,fill_value=(0.,1.),bounds_error=False)
        z_rnd = icdf(np.random.random(size=ntot))

        #Measure sigma_RSD from the catalog and apply it as a Gaussian shift.
        if DLA_catalog_path is None:
            raise ValueError('Needed a path to read the catalog')
        tab = fits.open(DLA_catalog_path)['DLACAT'].data
        z_master_NO_RSD = tab['Z_DLA_NO_RSD']
        z_master_RSD = tab['Z_DLA_RSD']
        dz_rsd_master = z_master_RSD - z_master_NO_RSD
        sigma_rsd_master = np.std(dz_rsd_master)
        z_rnd += np.random.normal(size=ntot,scale=sigma_rsd_master)

    dla_z = np.zeros(ntot)
    dla_skw_id = np.zeros(ntot,dtype='int32')
    dla_count = 0

    #For each DLA, place it in a skewer at random. Only keep it if it has
    #redshift lower than the quasar's.
    for i,dla_z_value in enumerate(z_rnd):
        skw_id = np.random.choice(n_qso)
        if dla_z_value < z_qso[skw_id]:
            dla_z[dla_count] = dla_z_value
            dla_skw_id[dla_count] = skw_id
            dla_count += 1
        print((i*100/ntot).round(5),end='\r')

    #Trim empty cells away.
    dla_z = dla_z[:dla_count]
    dla_skw_id = dla_skw_id[:dla_count]

    #Get the redshift of each DLA's skewer's QSO, its angular positions, MOCKID and pixel number.
    dla_ra = RA[dla_skw_id]
    dla_dec = DEC[dla_skw_id]
    dla_z_qso = z_qso[dla_skw_id]
    dla_z_qso_rsd = z_qso_rsd[dla_skw_id]
    dla_MOCKID = MOCKID[dla_skw_id]
    dla_pixnum = pixnum[dla_skw_id]

    #Make DLAIDs.
    dlaid = np.array(list(range(dla_count)))

    #Assign each DLA an NHI value if desired, and make a table.
    if add_NHI:
        dla_NHI = DLA.get_NHI(dla_z,NHI_min=NHI_min,NHI_max=NHI_max)
        dtype = [('RA', '>f8'), ('DEC', '>f8'), ('Z_QSO_NO_RSD', '>f8'), ('Z_QSO_RSD', '>f8'), ('Z_DLA', '>f8'), ('N_HI_DLA', '>f8'), ('MOCKID', '>i8'), ('DLAID', '>i8'), ('PIXNUM', '>i8')]
        DLA_data = np.array(list(zip(dla_ra,dla_dec,dla_z_qso,dla_z_qso_rsd,dla_z,dla_NHI,dla_MOCKID,dlaid,dla_pixnum)),dtype=dtype)
    else:
        dtype = [('RA', '>f8'), ('DEC', '>f8'), ('Z_QSO_NO_RSD', '>f8'), ('Z_QSO_RSD', '>f8'), ('Z_DLA', '>f8'), ('MOCKID', '>i8'), ('DLAID', '>i8'), ('PIXNUM', '>i8')]
        DLA_data = np.array(list(zip(dla_ra,dla_dec,dla_z_qso,dla_z_qso_rsd,dla_z,dla_MOCKID,dlaid,dla_pixnum)),dypte=dtype)

    #Write the file.
    DLA.write_DLA_master([DLA_data],out_path,N_side,overwrite=overwrite)

    return

# Execute
generate_rnd(factor=factor,out_path=out_path,DLA_catalog_path=DLA_catalog_path,QSO_catalog_path=QSO_catalog_path,footprint=footprint,lambda_min=lambda_min,lambda_max=lambda_max,NHI_min=NHI_min,NHI_max=NHI_max,overwrite=overwrite,N_side=N_side,add_NHI=add_NHI)
