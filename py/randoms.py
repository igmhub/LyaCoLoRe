import healpy as hp
import numpy as np
import warnings

from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import interp1d

from lyacolore import catalog, DLA, utils


def generate_rnd(factor=3, out_path= None, method='use_catalog', cat_path=None, footprint=None, nz_filename='input_files/Nz_qso_130618_2_colore1_hZs.txt', min_cat_z=1.8, max_cat_z=4.0, overwrite=False, N_side=16,start_MOCKID_rnd=10**10,seed=0):
    """
    Routine to generate a random catalog in 3D following
    certain N(z) distribution

    Args:
    ----
    factor: Size of the generated catalog (before masking)
    out_path: Name of output file where randoms will be saved (default: None)
    method: Method to generate the random catalog (default: 'random_choice')
    """

    #Set up a random sampler.
    state = np.random.RandomState(seed)

    #Set up vectors from n(z) file.
    N_vec = 500
    zedges = np.linspace(min_cat_z,max_cat_z,N_vec)
    zvec = (zedges[1:] + zedges[:-1])/2.
    dz = zvec[1] - zvec[0]

    #Generate random redshifts by one of 3 different methods:
    if method=='from_catalog':
        #Method 1: choose from the QSO catalog and add small deviations.
        zedges, dndz, ntot = get_dndz_from_catalog(cat_path,min_cat_z,max_cat_z,factor,hdu_name='CATALOG',z_name='Z_QSO_RSD')
        z_rnd = get_random_z(zedges,dndz,ntot,state,sig_noise=10**-6)

    elif method=='from_model':
        #Method 3: Calculate the cumulative distribution and interpolate.

        #Get the total number of QSOs.
        nz_file = Table.read(nz_filename,format='ascii')
        ntot = 4*180**2/np.pi*np.sum(nz_file['col2'])*(nz_file['col1'][1]-nz_file['col1'][0])
        ntot = int(factor*ntot)

        #Get dndz by interpolating input file.
        spl_z = interp1d(nz_file['col1'],nz_file['col2'],fill_value=0.)
        dndz = spl_z(zvec)

        # Get sigma_rsd to replicate effects of RSDs.
        sigma_rsd = get_sigma_dz_rsd(cat_path,hdu_name='CATALOG',z_norsd_name='Z_QSO_NO_RSD',z_rsd_name='Z_QSO_RSD')

        # Generate redshifts.
        z_rnd = get_random_z(zedges,dndz,ntot,state,sig_noise=sigma_rsd)

    #Assign random positions on the sky to the QSOs.
    #We have to iterate to ensure that the footprint is ok and that we have
    #enough QSOs.
    ra_rnd = 360.*state.uniform(size=len(z_rnd))
    cth_rnd = -1+2.*state.uniform(size=len(z_rnd))
    dec_rnd = np.arcsin(cth_rnd)*180/np.pi

    QSO_filter = utils.make_QSO_filter(footprint)
    good = QSO_filter(ra_rnd,dec_rnd)
    ra_rnd = ra_rnd[good]
    dec_rnd = dec_rnd[good]
    nrnd = ra_rnd.shape[0]
    while nrnd<ntot:
        extra_ra_rnd = 360.*state.uniform(size=ntot)
        extra_cth_rnd = -1+2.*state.uniform(size=ntot)
        extra_dec_rnd = np.arcsin(extra_cth_rnd)*180/np.pi
        good = QSO_filter(extra_ra_rnd,extra_dec_rnd)
        extra_ra_rnd = extra_ra_rnd[good]
        extra_dec_rnd = extra_dec_rnd[good]
        ra_rnd = np.concatenate((ra_rnd,extra_ra_rnd))
        dec_rnd = np.concatenate((dec_rnd,extra_dec_rnd))
        nrnd = ra_rnd.shape[0]
    if nrnd>ntot:
        ra_rnd = ra_rnd[:ntot]
        dec_rnd = dec_rnd[:ntot]

    #Assign MOCKIDs to the QSOs, checking that there's no overlap with those in
    #the master file.
    if cat_path is None:
        raise ValueError('Needed a path to read the catalog')
    tab = fits.open(cat_path)['CATALOG'].data
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

def generate_rnd_dla(factor=3, out_path=None , DLA_cat_path=None, QSO_cat_path=None, footprint='desi_pixel_plus', lambda_min=3470., lambda_max=6550., NHI_min=17.2, NHI_max=22.5, overwrite=False, N_side=16, add_NHI=True, method='from_catalog', start_DLAID_rnd=10**12, seed=0):
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
    zmin = lambda_min/utils.lya_rest - 1
    zmax = lambda_max/utils.lya_rest - 1
    N_vec = 500
    zedges = np.linspace(zmin,zmax,N_vec)
    zvec = (zedges[1:] + zedges[:-1])/2.
    dz = zvec[1] - zvec[0]

    #Get data about the QSO sample.
    h = fits.open(QSO_cat_path)
    QSO_data = h['CATALOG'].data
    QSO_data.sort(order='Z_QSO_NO_RSD')
    RA = QSO_data['RA']
    DEC = QSO_data['DEC']
    z_qso = QSO_data['Z_QSO_NO_RSD']
    z_qso_rsd = QSO_data['Z_QSO_RSD']
    pixnum = QSO_data['PIXNUM']
    MOCKID = QSO_data['MOCKID']
    n_qso = len(QSO_data)
    h.close()

    #Generate random redshifts for the DLAs.
    if method=='from_catalog':
        #Method 1: choose from the DLA catalog and add small deviations.
        zedges, dndz, ntot = get_dndz_from_catalog(DLA_cat_path,zmin,zmax,factor,hdu_name='DLACAT',z_name='Z_DLA_RSD')
        z_rnd = get_random_z(zedges,dndz,ntot,state,sig_noise=10**-6)

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

    elif method=='from_model':
        #Method 2: Calculate the cumulative distribution and interpolate.

        #Get dndz, and use it to get the mean number of DLAs in each cell.
        #Scale it up if necessary and calculate the total number of DLAs.
        dndz = DLA.dndz(zvec,NHI_min=NHI_min,NHI_max=NHI_max)
        mean_n = dz * dndz
        mean_n *= factor
        ntot = (np.sum(mean_n) * n_qso).astype('int')

        # Get sigma_rsd to replicate effects of RSDs.
        sigma_rsd = get_sigma_dz_rsd(DLA_cat_path,hdu_name='DLACAT',z_norsd_name='Z_DLA_NO_RSD',z_rsd_name='Z_DLA_RSD')

        # Generate redshifts.
        z_rnd = get_random_z(zedges,dndz,ntot,state,sig_noise=sigma_rsd)

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
    if DLA_cat_path is None:
        raise ValueError('Needed a path to read the catalog')
    tab = fits.open(DLA_cat_path)['DLACAT'].data
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

def get_dndz_from_catalog(cat_path,zmin,zmax,factor,nvec=500,pc_extra=0.05,hdu_name='CATALOG',z_name='Z_QSO_NO_RSD'):

    #Make a z vector with edges extended from min and max cat values.
    #This is to take into account RSDs in QSOs near the edge of the range.
    extra = (zmax - zmin) * pc_extra
    zedges = np.linspace(zmin-extra,zmax+extra,nvec*(1+2*extra))
    zvec = (zedges[1:] + zedges[:-1])/2.
    dz = zvec[1] - zvec[0]

    #Get catalog data and calculate dndz from it.
    if cat_path is None:
        raise ValueError('Need a path to read the catalog.')
    tab = fits.open(cat_path)[hdu_name].data
    z_cat = tab[z_name]
    n_cat = len(z_cat)
    dndz, _ = np.histogram(z_cat,bins=zedges)

    #Get ntot
    ntot = int(n_cat * factor)

    return zedges, dndz, ntot

def get_random_z(zedges,dndz,ntot,state,sig_noise=10**-6):

    cdf = np.cumsum(dndz)/np.sum(dndz)
    cdf_i = np.concatenate([[0],cdf])
    icdf = interp1d(cdf_i,zedges,fill_value=(0.,1.),bounds_error=False)
    z_rnd = icdf(state.uniform(size=ntot)) + state.normal(size=ntot,scale=sig_noise)

    return z_rnd

def get_sigma_dz_rsd(cat_path,hdu_name='CATALOG',z_norsd_name='Z_QSO_NO_RSD',z_rsd_name='Z_QSO_RSD'):

    tab = fits.open(cat_path)[hdu_name].data
    z_norsd = tab[z_norsd_name]
    z_rsd = tab[z_rsd_name]
    dz_rsd = z_rsd - z_norsd
    sigma_dz_rsd = np.std(dz_rsd)

    return sigma_dz_rsd
