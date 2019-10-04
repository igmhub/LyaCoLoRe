import numpy as np
from scipy.interpolate import interp1d
import astropy.io.fits as fits
import warnings

from lyacolore import DLA, utils

lya = utils.lya_rest

#Set up options
factor = 10.
#basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v9/v9.0.0_full/'
#basedir = '/global/projecta/projectdirs/desi/mocks/lya_forest/develop/london/v8.0/v8.0.0/'
basedir = '/project/projectdirs/desi/users/jfarr/LyaCoLoRe_paper/data/LyaCoLoRe_output/v9.0.9/'
out_path = basedir+'/master_DLA_randoms.fits'
method = 'from_catalog'
DLA_catalog_path = basedir+'/master_DLA.fits'
QSO_catalog_path = basedir+'/master.fits'
footprint = 'desi_pixel_plus'
lambda_min = 3470.
lambda_max = 6550.
NHI_min = 17.2
NHI_max = 22.5
overwrite = True
N_side = 16
add_NHI = False
start_DLAID_rnd = 10**12

def generate_rnd(factor=3, out_path=None , DLA_catalog_path=None, QSO_catalog_path=None, footprint='desi_pixel_plus', lambda_min=3470., lambda_max=6550., NHI_min=17.2, NHI_max=22.5, overwrite=False, N_side=16, add_NHI=True, method='from_catalog', start_DLAID_rnd=10**12, seed=0):
    """
    Routine to generate a random catalog in 3D following
    certain N(z) distribution

    Args:
    ----
    factor: Size of the generated catalog (before masking)
    out_path: Output path
    """

    #Set up a random sampler.
    state = np.random.RandomState(seed)

    #Generate a z vector and the dn/dz function.
    zmin = lambda_min/lya - 1
    zmax = lambda_max/lya - 1
    N_vec = 500
    zvec = np.linspace(zmin,zmax,N_vec)
    zedges = np.concatenate([[zvec[0]]-(zvec[1]-zvec[0])/2.,(zvec[1:]+zvec[:-1])*0.5,[zvec[-1]+(-zvec[-2]+zvec[-1])*0.5]]).ravel()
    dz = zvec[1] - zvec[0]

    #Get data about the QSO sample.
    h = fits.open(QSO_catalog_path)
    QSO_data = h['CATALOG'].data
    QSO_data.sort(order='Z_QSO_NO_RSD')
    RA = QSO_data['RA']
    DEC = QSO_data['DEC']
    z_qso = QSO_data['Z_QSO_NO_RSD']
    z_qso_rsd = QSO_data['Z_QSO_RSD']
    pixnum = QSO_data['PIXNUM']
    MOCKID = QSO_data['MOCKID']
    n_qso = z_qso.shape[0]
    h.close()

    #Generate random redshifts for the DLAs.
    if method=='from_catalog':
        #Method 1: choose from the DLA catalog and add small deviations.
        #Make a z vector with edges extended from min and max cat values.
        #This is to take into account RSDs in QSOs near the edge of the range.
        pc_extra = 0.05
        extra = (zmax - zmin) * pc_extra
        zedges = np.linspace(zmin-extra,zmax+extra,N_vec*(1+2*extra))
        zvec = (zedges[1:] + zedges[:-1])/2.
        dz = zvec[1] - zvec[0]

        #Get catalog data and calculate dndz from it.
        if DLA_catalog_path is None:
            raise ValueError('Needed a path to read the catalog')
        tab = fits.open(DLA_catalog_path)['DLACAT'].data
        z_master_RSD = tab['Z_DLA_RSD']
        n_DLA_cat = z_master_RSD.shape[0]
        dndz_RSD,_ = np.histogram(z_master_RSD,bins=zedges)

        #Calculate n_total from the number in the catalog.
        ntot = int(n_DLA_cat * factor)

        #Turn dn/dz into a cdf and draw redshifts from it.
        cdf_RSD = np.cumsum(dndz_RSD)/np.sum(dndz_RSD)
        cdf_RSD_i = np.concatenate([[0],cdf_RSD])
        icdf_RSD = interp1d(cdf_RSD_i,zedges,fill_value=(0.,1.),bounds_error=False)
        z_rnd = icdf_RSD(state.uniform(size=ntot))

        #Exclude DLAs that are beyond the furthest QSO.
        #In reality should these be excluded?
        w = z_rnd<np.max(z_qso)
        z_rnd = z_rnd[w]
        dla_z = np.zeros(ntot)
        dla_skw_id = np.zeros(ntot,dtype='int32')
        dla_count = 0

        #For each DLA, place it in a skewer at random.
        #As we are using a *measured* dn/dz, we cannot discard any. Thus, for
        #each DLA we search for a valid QSO: i.e. one with z_QSO > z_DLA.
        max_n_att = 100
        for i,dla_z_value in enumerate(z_rnd):
            valid = False
            n_att = 0
            low = 0
            while valid == False:
                #Choose a random skewer.
                skw_id = state.randint(low=low,high=n_qso)
                #print(low,n_qso,skw_id)
                n_att += 1
                #If the QSO is valid, stop searching.
                if z_qso[skw_id] > dla_z_value:
                    valid = True
                #Otherwise, if we have exceeded the maximum number of random
                #attempts, find which skewers are valid.
                elif n_att > max_n_att:
                    print('\nMAX N ATTEMPTS REACHED\n')
                    loc = np.searchsorted(z_qso,dla_z_value)
                    #If there exist valid skewers, choose one of them.
                    if loc < n_qso:
                        skw_id = state.randint(low=loc,high=n_qso)
                        valid = True
                    #Otherwise exit.
                    else:
                        break
                else:
                    #print('\nSkewer',skw_id,'is too low at',z_qso[skw_id],'for dla with z',dla_z_value,'. Look at',skw_id+1,'out of',n_qso)
                    low = skw_id+1
            #If a valid skewer was found, store the information.
            if valid:
                dla_z[dla_count] = dla_z_value
                dla_skw_id[dla_count] = skw_id
                dla_count += 1
            print(round(i*100/ntot,5),end='\r')

        #Trim empty cells away.
        dla_z = dla_z[:dla_count]
        dla_skw_id = dla_skw_id[:dla_count]

    elif method=='cdf':
        #Method 2: Calculate the cumulative distribution and interpolate.
        #Get dndz, and use it to get the mean number of DLAs in each cell.
        #Scale it up if necessary and calculate n_total.
        dndz = DLA.dndz(zvec,NHI_min=NHI_min,NHI_max=NHI_max)
        mean_n = dz * dndz
        mean_n *= factor
        ntot = (np.sum(mean_n) * n_qso).astype('int')
        print(ntot)

        #Generate redshifts without RSDs.
        cdf = np.cumsum(dndz)/np.sum(dndz)
        cdf_i = np.concatenate([[0],cdf])
        icdf = interp1d(cdf_i,zedges,fill_value=(0.,1.),bounds_error=False)
        z_rnd = icdf(state.uniform(size=ntot))

        #Measure sigma_RSD from the catalog and apply it as a Gaussian shift.
        if DLA_catalog_path is None:
            raise ValueError('Needed a path to read the catalog')
        tab = fits.open(DLA_catalog_path)['DLACAT'].data
        z_master_NO_RSD = tab['Z_DLA_NO_RSD']
        z_master_RSD = tab['Z_DLA_RSD']
        dz_rsd_master = z_master_RSD - z_master_NO_RSD
        sigma_rsd_master = np.std(dz_rsd_master)
        z_rnd += state.normal(size=ntot,scale=sigma_rsd_master)

        dla_z = np.zeros(ntot)
        dla_skw_id = np.zeros(ntot,dtype='int32')
        dla_count = 0

        #For each DLA, place it in a skewer at random. Only keep it if it has
        #redshift lower than the QSO's.
        for i,dla_z_value in enumerate(z_rnd):
            skw_id = state.choice(n_qso)
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

    #Assign DLAIDs to the DLAs, checking that there's no overlap with those in
    #the master file.
    if DLA_catalog_path is None:
        raise ValueError('Needed a path to read the catalog')
    tab = fits.open(DLA_catalog_path)['DLACAT'].data
    max_cat_DLAID = np.max(tab['DLAID'])
    while max_cat_DLAID > start_DLAID_rnd:
        warnings.warn('Start value of randoms\' MOCKIDs is not high enough: increasing from {} to {}'.format(start_DLAID_rnd,10*start_DLAID_rnd))
        start_DLAID_rnd *= 10
    dlaid = (np.array(list(range(dla_count))) + start_DLAID_rnd).astype('int')
    if np.max(dlaid) > (2**63 - 1):
        raise ValueError('Max DLAID exceeds max integer allowed by FITS.')

    #Assign each DLA an NHI value if desired, and make a table.
    if add_NHI:
        print('\nGetting NHI values...')
        dla_NHI = DLA.get_NHI(dla_z,NHI_min=NHI_min,NHI_max=NHI_max)
        dtype = [('RA', '>f8'), ('DEC', '>f8'), ('Z_QSO_NO_RSD', '>f8'), ('Z_QSO_RSD', '>f8'), ('Z_DLA', '>f8'), ('N_HI_DLA', '>f8'), ('MOCKID', '>i8'), ('DLAID', '>i8'), ('PIXNUM', '>i8')]
        DLA_data = np.array(list(zip(dla_ra,dla_dec,dla_z_qso,dla_z_qso_rsd,dla_z,dla_NHI,dla_MOCKID,dlaid,dla_pixnum)),dtype=dtype)
    else:
        dtype = [('RA', '>f8'), ('DEC', '>f8'), ('Z_QSO_NO_RSD', '>f8'), ('Z_QSO_RSD', '>f8'), ('Z_DLA', '>f8'), ('MOCKID', '>i8'), ('DLAID', '>i8'), ('PIXNUM', '>i8')]
        DLA_data = np.array(list(zip(dla_ra,dla_dec,dla_z_qso,dla_z_qso_rsd,dla_z,dla_MOCKID,dlaid,dla_pixnum)),dtype=dtype)

    #Write the file.
    DLA.write_DLA_master([DLA_data],out_path,N_side,overwrite=overwrite)

    return

# Execute
generate_rnd(factor=factor,out_path=out_path,DLA_catalog_path=DLA_catalog_path,QSO_catalog_path=QSO_catalog_path,footprint=footprint,lambda_min=lambda_min,lambda_max=lambda_max,NHI_min=NHI_min,NHI_max=NHI_max,overwrite=overwrite,N_side=N_side,add_NHI=add_NHI,method=method,start_DLAID_rnd=start_DLAID_rnd,seed=i)
