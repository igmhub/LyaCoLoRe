from __future__ import print_function, division
import numpy as np
import astropy.io.fits as fits
from scipy.stats import norm
from scipy.interpolate import interp1d, interp2d
import astropy.table
import os
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15
import time

from lyacolore import utils

try:
    ## Use Quiet class to stop pyigm from printing output.
    with utils.Quiet():
        from pyigm.fN.fnmodel import FNModel
        fN_default = FNModel.default_model(cosmo=Planck15)
    fN_default.zmnx = (0.,5.)
    fN_cosmo = fN_default.cosmo
    use_pyigm = True
except:
    print('WARN: pyigm not found, DLA distributions will not be perfect!')
    use_pyigm = False

def get_threshold(b2f):
    """ Compute the field threshold for a given bias
    Arguments:
     - b2f (array): the bias of DLAs relative to the field of interest

    Outputs:
     - threshold (array): for each input bias value, the threshold value of
     the field
    """
    nu = np.linspace(-10,100,500) # Generous range to interpolate
    p_nu = norm.pdf(nu)
    galaxy_mean = 1.0-norm.cdf(nu)
    bDs_nu = np.zeros(nu.shape)
    bDs_nu[galaxy_mean!=0] = p_nu[galaxy_mean!=0]/galaxy_mean[galaxy_mean!=0]
    y = interp1d(bDs_nu,nu)
    threshold = y(b2f)
    return threshold

def get_bias_z(fname,dla_bias):
    """ Given a path, read the z array there and return a bias inversely
    proportional to the growth"""
    colore_cosmo = fits.open(fname)[4].data
    z = colore_cosmo['Z']
    D = colore_cosmo['D']
    y = interp1d(z,D)
    bias = dla_bias/D*y(2.25)
    return z, bias, D

def get_sigma_g(object,mode='global'):
    if mode=='global':
        #Use the 'global' sigma_G - that which is not biased by cells close to
        #the QSOs. This is constant with redshift.
        return object.SIGMA_G
    if mode=='sample':
        #Use the average sigma_G with redshift, including cells near to QSOs.
        return object.sample_SIGMA_G
    else:
        raise ValueError('DLA bias method not recognised.')

def flag_DLA(z_qso,z_cells,deltas,nu_arr,sigma_g):
    """ Flag the pixels in a skewer where DLAs are possible"""
    # find cells with density above threshold
    # TODO: why do we multiply nu by sigma_g here?
    flag = deltas > nu_arr*sigma_g
    # mask cells with z > z_qso, where DLAs would not be observed
    Nq=len(z_qso)
    for i in range(Nq):
        low_z = z_cells < z_qso[i]
        flag[i,:] *= low_z
    return flag

#number per unit redshift from minimum lg(N) in file (17.2) to argument
# Reading file from https://arxiv.org/pdf/astro-ph/0407378.pdf

def dnHD_dz_cumlgN(z,logN):
    tab = astropy.table.Table.read(os.environ['LYACOLORE_PATH']+'/input_files/dla_files/zheng_cumulative.overz',format='ascii')
    y = interp2d(tab['col1'],tab['col2'],tab['col3'],fill_value=None)
    return y(z,logN)

def dndz(z, NHI_min=17.2, NHI_max=22.5):
    """ Get the column density distribution as a function of z,
        for a given range in N"""

    if use_pyigm:
        # get incidence rate per path length dX (in comoving coordinates)
        dndX = fN_default.calculate_lox(z,NHI_min,NHI_max)
        # convert dX to dz
        dXdz = fN_cosmo.abs_distance_integrand(z)
        dndz = dndX * dXdz
    else:
        dndz = dnHD_dz_cumlgN(z,NHI_max)-dnHD_dz_cumlgN(z,NHI_min)
 
    return dndz

def get_NHI(z, NHI_min=17.2, NHI_max=22.5, NHI_nsamp=100):
    """ Get random column densities for a given z
    """
    times = []
    t = time.time()

    # number of DLAs we want to generate
    Nz = len(z)

    # Set up the grid in NHI, and define its edges/widths.
    # First in log space.
    log_NHI_edges = np.linspace(NHI_min,NHI_max,NHI_nsamp+1)
    log_NHI = (log_NHI_edges[1:] + log_NHI_edges[:-1])/2.
    log_NHI_widths = log_NHI_edges[1:] - log_NHI_edges[:-1]
    # Then in linear space.
    NHI_edges = 10**log_NHI_edges
    NHI_widths = NHI_edges[1:] - NHI_edges[:-1]

    times += [time.time()-t]
    t = time.time()

    probs = np.zeros([Nz,NHI_nsamp])

    if use_pyigm:

        #Evaluate f at the points of the NHI grid and each redshift.
        f = 10**fN_default.evaluate(log_NHI,z)

        times += [time.time()-t]
        t = time.time()

        #Calculate the probaility of each NHI bin.
        aux = f*np.outer(NHI_widths,np.ones(z.shape))

        times += [time.time()-t]
        t = time.time()

        probs = (aux/np.sum(aux,axis=0)).T

        times += [time.time()-t]
        t = time.time()

    else:

        # TODO: test this
        probs_low = dnHD_dz_cumlgN(z,log_NHI_edges[:-1]).T
        probs_high = dnHD_dz_cumlgN(z,log_NHI_edges[1:]).T
        probs = probs_high-probs_low

    times += [time.time()-t]
    t = time.time()

    #Calculate the cumulative distribution
    """
    cumulative = np.zeros(probs.shape)
    for i in range(NHI_nsamp):
        cumulative[:,i] = np.sum(probs[:,:i],axis=1)
    """
    cumulative = np.cumsum(probs,axis=1)

    #Add the top and bottom edges on to improve interpolation.
    """
    log_NHI_interp = np.concatenate([[log_NHI_edges[0]],log_NHI,[log_NHI_edges[1]]])
    end_0 = np.zeros((z.shape[0],1))
    end_1 = np.ones((z.shape[0],1))
    cumulative_interp = np.concatenate([end_0,cumulative,end_1],axis=1)
    """

    log_NHI_interp = log_NHI_edges
    end_0 = np.zeros((z.shape[0],1))
    cumulative_interp = np.concatenate([end_0,cumulative],axis=1)

    times += [time.time()-t]
    t = time.time()

    #Assign NHI values by choosing a random number in [0,1] and interpolating
    #the cumulative distribution to get a value of NHI.
    log_NHI_values = np.zeros(Nz)
    for i in range(Nz):
        p = np.random.uniform()
        log_NHI_values[i] = np.interp(p,cumulative_interp[i,:],log_NHI_interp)

    times += [time.time()-t]
    t = time.time()

    times = np.array(times)
    times /= np.sum(times)
    #print(times.round(4))

    return log_NHI_values

def get_DLA_table(object,dla_bias=2.0,dla_bias_z=2.25,extrapolate_z_down=None,NHI_min=17.2,NHI_max=22.5,seed=123,evol='b_const',method='global'):

    #Hopefully this sets the seed for all random generators used
    np.random.seed(seed)

    #Extract data from the object.
    z_qso = object.Z_QSO
    z_cell = object.Z
    D_cell = object.D

    #Get the edges of the z cells, and their widths
    if extrapolate_z_down and extrapolate_z_down<z_cell[0]:
        zedges = np.concatenate([[extrapolate_z_down],(z_cell[1:]+z_cell[:-1])*0.5,[z_cell[-1]+(-z_cell[-2]+z_cell[-1])*0.5]]).ravel()
    else:
        zedges = np.concatenate([[z_cell[0]]-(z_cell[1]-z_cell[0])/2.,(z_cell[1:]+z_cell[:-1])*0.5,[z_cell[-1]+(-z_cell[-2]+z_cell[-1])*0.5]]).ravel()
    z_width = zedges[1:]-zedges[:-1]

    #Setup bias as a function of redshift: either b constant with z, or b*D constant with z.
    y = interp1d(z_cell,D_cell)

    ## Determine the threshold above which to flag cells as potential DLA
    ## locations. This depends on the skewer type and the DLA bias evolution.
    if object.GAUSSIAN_DELTA_rows is not None:
        sigma_g = get_sigma_g(object,mode=method)
        if evol == "b_const":
            b2f = dla_bias*D_cell*sigma_g
        elif evol == "bD_const":
            b2f = dla_bias*y(dla_bias_z)*sigma_g*np.ones(z.shape)
        else:
            raise ValueError('DLA bias evol not recognised.')
        deltas = object.GAUSSIAN_DELTA_rows

    elif object.DENSITY_DELTA_rows is not None:
        """
        ## Also need the pdf pf the skewers.
        ## For now, just raise an error.
        if evol == "b_const":
            b2f = dla_bias
        elif evol == "bD_const":
            b2f = dla_bias*y(dla_bias_z)*np.ones(z.shape)/D_cell
        else:
            raise ValueError('DLA bias evol not recognised.')
        deltas = object.DENSITY_DELTA_rows
        """
        raise ValueError('Adding DLAs without Gaussian skewers is not implemented as yet.')

    else:
        raise ValueError('Cannot find Gaussian or density delta skewers when flagging cells for DLAs')

    nu_arr = get_threshold(b2f)
    flagged_cells = flag_DLA(z_qso,z_cell,deltas,nu_arr,sigma_g)

    #Mean of the poisson is the mean number of DLAs with redshift scaled up by
    #the proportion of expected flagged cells.
    mean_N_per_cell = z_width * dndz(z_cell,NHI_min=NHI_min,NHI_max=NHI_max)
    p_nu_z = 1.0-norm.cdf(nu_arr)
    w = p_nu_z>0
    mu = np.zeros(p_nu_z.shape)
    mu[w] = mean_N_per_cell[w]/p_nu_z[w]

    #Draw number of DLAs per cell from a Poisson distribution, place them only
    #in flagged cells.
    pois = np.random.poisson(mu,size=(len(z_qso),len(mu)))
    dlas_in_cell = pois*flagged_cells

    #Store information for each of the DLAs that will be added
    ndlas = np.sum(dlas_in_cell)
    dla_z = np.zeros(ndlas)
    dla_skw_id = np.zeros(ndlas, dtype='int32')
    dla_rsd_dz = np.zeros(ndlas)
    dla_count = 0

    #In each skewer, identify cells where we need to add a DLA
    for skw_id,dla in enumerate(dlas_in_cell):

        #Find cells that will be allocated at least one DLA
        dla_cells = np.where(dla>0)[0]

        #For each dla, assign it a redshift, a velocity and a column density.
        for cell in dla_cells:
            dla_z[dla_count:dla_count+dla[cell]] = np.random.uniform(low=(zedges[cell]),high=(zedges[cell+1]),size=dla[cell])
            dla_skw_id[dla_count:dla_count+dla[cell]] = skw_id
            dla_rsd_dz[dla_count:dla_count+dla[cell]] = object.VEL_rows[skw_id,cell]
            dla_count = dla_count+dla[cell]

    dla_NHI = get_NHI(dla_z,NHI_min=NHI_min,NHI_max=NHI_max)

    #Obtain the other variable to put in the DLA table
    MOCKIDs = object.MOCKID[dla_skw_id]

    #Make DLAIDs
    DLAIDs = np.zeros(MOCKIDs.shape)
    current_MOCKID = 0
    current_DLAID = current_MOCKID * 10**3
    for i,MOCKID in enumerate(MOCKIDs):
        if MOCKID != current_MOCKID:
            current_MOCKID = MOCKID
            current_DLAID = current_MOCKID * 10**3
        DLAIDs[i] = current_DLAID.astype('int')
        current_DLAID += 1

    #Make the data into a table HDU
    data = list(zip(dla_z,dla_z+dla_rsd_dz,dla_NHI,MOCKIDs,DLAIDs))
    dtype = [('Z_DLA_NO_RSD', '>f8'), ('Z_DLA_RSD', '>f8'), ('N_HI_DLA', '>f8'), ('MOCKID', '>i8'), ('DLAID', '>i8')]
    dla_table = np.array(data,dtype=dtype)

    return dla_table

def get_DLA_data_from_transmission(pixel,filename):

    DLA_data = []
    t = fits.open(filename)
    DLA_table = np.sort(t['DLA'].data,order=['DLAID'])

    current_MOCKID = 0
    #current_DLAID = current_MOCKID * 10**3

    for i,DLA in enumerate(DLA_table):
        MOCKID = DLA['MOCKID']
        DLAID = DLA['DLAID']

        Z_DLA_NO_RSD = DLA['Z_DLA_NO_RSD']
        try:
            Z_DLA_RSD = DLA['Z_DLA_RSD']
        except KeyError:
            #This is to allow for a typo made in v4.2, that labelled Z_DLA wrongly. It should be removed afterwards.
            Z_DLA_RSD = DLA['DZ_DLA_RSD']
        N_HI_DLA = DLA['N_HI_DLA']

        if MOCKID != current_MOCKID:

            QSO_data = t[1].data[t[1].data['MOCKID']==MOCKID]
            RA = QSO_data['RA']
            DEC = QSO_data['DEC']
            Z_QSO_RSD = QSO_data['Z']
            Z_QSO_NO_RSD = QSO_data['Z_noRSD']

            current_MOCKID = MOCKID
            #current_DLAID = current_MOCKID * 10**3

        #DLAID = current_DLAID
        #current_DLAID += 1

        DLA_data += [(RA,DEC,Z_QSO_NO_RSD,Z_QSO_RSD,Z_DLA_NO_RSD,Z_DLA_RSD,N_HI_DLA,MOCKID,DLAID,pixel)] #No file number

    t.close()

    dtype = [('RA', '>f8'), ('DEC', '>f8'), ('Z_QSO_NO_RSD', '>f8'), ('Z_QSO_RSD', '>f8'), ('Z_DLA_NO_RSD', '>f8'), ('Z_DLA_RSD', '>f8'), ('N_HI_DLA', '>f8'), ('MOCKID', '>i8'), ('DLAID', '>i8'), ('PIXNUM', '>i8')]
    DLA_data = np.array(DLA_data,dtype=dtype)

    return DLA_data

def write_DLA_master(DLA_data_list,filename,N_side,overwrite=False,picca_DRQ=False,add_DRQ_RSDs=True):

    DLA_master_data = np.concatenate(DLA_data_list)

    if picca_DRQ:
        if add_DRQ_RSDs:
            DRQ_RSD_option = 'RSD'
        else:
            DRQ_RSD_option = 'NO_RSD'
        RA = DLA_master_data['RA']
        DEC = DLA_master_data['DEC']
        Z = DLA_master_data['Z_DLA_'+DRQ_RSD_option]
        THING_ID = DLA_master_data['MOCKID']
        MJD = np.zeros(RA.shape[0])
        FID = np.zeros(RA.shape[0])
        PLATE = THING_ID
        PIXNUM = DLA_master_data['PIXNUM']

        #Make the data array.
        dtype = [('RA','f4'),('DEC','f4'),('Z','f4'),('THING_ID',int),('MJD','f4'),('FIBERID',int),('PLATE',int),('PIXNUM',int)]
        DLA_master_data = np.array(list(zip(RA,DEC,Z,THING_ID,MJD,FID,PLATE,PIXNUM)),dtype=dtype)

    #Make an appropriate header.
    header = fits.Header()
    header['NSIDE'] = N_side

    #Make the data into tables.
    hdu_ID = fits.BinTableHDU.from_columns(DLA_master_data,header=header,name='DLACAT')

    #Make a primary HDU.
    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)

    #Make the .fits file.
    hdulist = fits.HDUList([prihdu,hdu_ID])
    hdulist.writeto(filename,overwrite=overwrite)
    hdulist.close()

    return
