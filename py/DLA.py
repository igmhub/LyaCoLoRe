from __future__ import print_function, division
import numpy as np
import astropy.io.fits as fits
from scipy.stats import norm
from scipy.interpolate import interp1d, interp2d
import astropy.table
import os
import matplotlib.pyplot as plt

try:
    from pyigm.fN.fnmodel import FNModel
    fN_default = FNModel.default_model()
    fN_default.zmnx = (0.5,4)
    fN_cosmo = fN_default.cosmo
    use_pyigm = True
except:
    use_pyigm = False

def nu_of_bD(b):
    """ Compute the Gaussian field threshold for a given bias"""
    nu = np.linspace(-10,100,500) # Generous range to interpolate
    p_nu = norm.pdf(nu)
    galaxy_mean = 1.0-norm.cdf(nu)
    b_nu = np.zeros(nu.shape)
    b_nu[galaxy_mean!=0] = p_nu[galaxy_mean!=0]/galaxy_mean[galaxy_mean!=0]
    y = interp1d(b_nu,nu)
    return y(b)

def get_bias_z(fname,dla_bias):
    """ Given a path, read the z array there and return a bias inversely
    proportional to the growth"""
    colore_cosmo = fits.open(fname)[4].data
    z = colore_cosmo['Z']
    D = colore_cosmo['D']
    y = interp1d(z,D)
    bias = dla_bias/D*y(2.25)
    return z, bias, D

def get_sigma_g(fname, mode='SG'):
    if mode=='SG':
        # Biased as well
        return fits.open(fname)[4].header['SIGMA_G']
    if mode=='SKW':
        # Approximation 2: Take the skewers (biased when QSOs are present)
        skewers = fits.open(fname)[2].data
        return np.std(skewers,axis=0)

def flag_DLA(zq,z_cells,deltas,nu_arr,sigma_g):
    """ Flag the pixels in a skewer where DLAs are possible"""
    # find cells with density above threshold
    flag = deltas > nu_arr*sigma_g
    # mask cells with z > z_qso, where DLAs would not be observed
    Nq=len(zq)
    for i in range(Nq):
        low_z = z_cells < zq[i]
        flag[i,:] *= low_z
    return flag

#number per unit redshift from minimum lg(N) in file (17.2) to argument
# Reading file from https://arxiv.org/pdf/astro-ph/0407378.pdf

def dnHD_dz_cumlgN(z,logN):
    tab = astropy.table.Table.read(os.path.abspath('example_data/zheng_cumulative.overz'),format='ascii')
    y = interp2d(tab['col1'],tab['col2'],tab['col3'],fill_value=None)
    return y(z,logN)

def dNdz(z, Nmin=19.5, Nmax=22.):
    """ Get the column density distribution as a function of z,
    for a given range in N"""
    if use_pyigm:
        # get incidence rate per path length dX (in comoving coordinates)
        dNdX = fN_default.calculate_lox(z,Nmin,Nmax)
        # convert dX to dz
        dXdz = fN_cosmo.abs_distance_integrand(z)
        return dNdX * dXdz
    else:
        return dnHD_dz_cumlgN(z,Nmax)-dnHD_dz_cumlgN(z,Nmin)

def get_N(z, Nmin=19.5, Nmax=22.0, nsamp=100):
    """ Get random column densities for a given z
    """
   
    # number of DLAs we want to generate
    Nz = len(z) 
    nn = np.linspace(Nmin,Nmax,nsamp)
    probs = np.zeros([Nz,nsamp])
    if use_pyigm:
        auxfN = np.cumsum(fN_default.evaluate(nn,z), axis=0).T
        probs_low = auxfN[:,1:]
        probs_high = auxfN[:,:-1]
    else:
        probs_low = dnHD_dz_cumlgN(z,nn[:-1]).T 
        probs_high = dnHD_dz_cumlgN(z,nn[1:]).T 
    probs[:,1:] = probs_high-probs_low
    NHI = np.zeros(Nz)
    for i in range(Nz):
        NHI[i] = np.random.choice(nn,size=1,p=probs[i]/np.sum(probs[i]))+(nn[1]-nn[0])*np.random.random(size=1)
    return NHI

def add_DLA_table_to_object(object,dla_bias=2.0,extrapolate_z_down=None,Nmin=19.5,Nmax=22.,seed=123):

    #Hopefully this sets the seed for all random generators used
    np.random.seed(seed)

    #Quasar redshift for each skewer
    zq = object.Z_QSO
    #Redshift of each cell in the skewer
    z_cell = object.Z
    #Linear growth rate of each cell in the skewer
    D_cell = object.D

    #Setup bias as a function of redshift
    y = interp1d(z_cell,D_cell)
    bias = dla_bias/(D_cell)*y(2.25)

    #We measure sigma_G already, but it is not fed back into the files at all. This should change.
    sigma_g = object.SIGMA_G
    #sigma_g = DLA.get_sigma_g(o.input_file)

    nu_arr = nu_of_bD(bias*D_cell)
    #Figure out cells that could host a DLA, based on Gaussian fluctuation
    deltas = object.GAUSSIAN_DELTA_rows
    flagged_cells = flag_DLA(zq,z_cell,deltas,nu_arr,sigma_g)

    #Edges of the z bins
    if extrapolate_z_down and extrapolate_z_down<z_cell[0]:
        zedges = np.concatenate([[extrapolate_z_down],(z_cell[1:]+z_cell[:-1])*0.5,[z_cell[-1]+(-z_cell[-2]+z_cell[-1])*0.5]]).ravel()
    else:
        zedges = np.concatenate([[z_cell[0]],(z_cell[1:]+z_cell[:-1])*0.5,[z_cell[-1]+(-z_cell[-2]+z_cell[-1])*0.5]]).ravel()
    z_width = zedges[1:]-zedges[:-1]

    #Get the average number of DLAs per cell, from the column density dist.
    mean_N_per_cell = z_width*dNdz(z_cell,Nmin=Nmin,Nmax=Nmax)

    #For a given z, probability of having the density higher than the threshold
    p_nu_z = 1.0-norm.cdf(nu_arr)

    #Define mean of the Poisson distribution (per cell)
    mu = mean_N_per_cell/p_nu_z

    #Select cells that will hold a DLA, drawing from the Poisson distribution
    pois = np.random.poisson(mu,size=(len(zq),len(mu)))
    #Number of DLAs in each cell (mostly 0, several 1, not many with >1)
    dlas_in_cell = pois*flagged_cells

    #Total number of DLAs to be added to all the cells of all the skewers
    ndlas = np.sum(dlas_in_cell)
    #Store information for each of the DLAs that will be added
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

    dla_NHI = get_N(dla_z,Nmin=Nmin,Nmax=Nmax)



    #global id for the skewers
    MOCKIDs = object.MOCKID[dla_skw_id]

    #Make the data into a table HDU
    dla_table = astropy.table.Table([MOCKIDs,dla_z,dla_rsd_dz,dla_NHI],names=('MOCKID','Z_DLA','DZ_DLA','N_HI_DLA'))

    #print('DLA table',dla_table)

    object.DLA_table = dla_table

