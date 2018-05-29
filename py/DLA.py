from __future__ import print_function, division
import numpy as np
import astropy.io.fits as fits
from scipy.stats import norm
from scipy.interpolate import interp1d, interp2d
import astropy.table
import os

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
    cosmo = fits.open(fname)[4].data
    z = cosmo['Z']
    D = cosmo['D']
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

def flag_DLA(skewer,nu_arr,sigma_g):
    """ Flag the pixels in a skewer where DLAs are possible"""
    flag = skewer > nu_arr*sigma_g
    return flag

#number per unit redshift from minimum lg(N) in file (17.2) to argument
# Reading file from https://arxiv.org/pdf/astro-ph/0407378.pdf

def dnHD_dz_cumlgN(z,logN):
    tab = astropy.table.Table.read(os.path.abspath('example_data/zheng_cumulative.overz'),format='ascii')
    y = interp2d(tab['col1'],tab['col2'],tab['col3'],fill_value=None)
    return y(z,logN)

def dNdz(z, Nmin=19.5, Nmax=22.):
    """ Get the column density as a function of z
    for a given range in N"""
    return dnHD_dz_cumlgN(z,Nmax)-dnHD_dz_cumlgN(z,Nmin)

def get_N(z, Nmin=19.5, Nmax=22.0, nsamp=100):
    """ Get the column density for a given z
    This always returns recurring decimals of a kind, could just expand nsamp to deal with it"""
    nn = np.linspace(Nmin,Nmax,nsamp)
    probs = dnHD_dz_cumlgN(z,nn).T
    N = np.zeros(len(probs))
    for i in range(0,len(probs)):
        N[i] = np.random.choice(nn,size=1,p=probs[i]/np.sum(probs[i]))
    return N

def add_DLA_table_to_object(object,dla_bias=2.0):

    y = interp1d(object.Z,object.D)
    bias = dla_bias/(object.D)*y(2.25)

    #We measure sigma_G already, but it is not fed back into the files at all. This should change.
    sigma_g = object.SIGMA_G
    #sigma_g = DLA.get_sigma_g(o.input_file)

    nu_arr = nu_of_bD(bias*object.D)
    flagged_pixels = flag_DLA(object.GAUSSIAN_DELTA_rows,nu_arr,sigma_g)

    #Edges of the z bins
    zedges = np.concatenate([[object.Z[0]],(object.Z[1:]+object.Z[:-1])*0.5,[object.Z[-1]+(-object.Z[-2]+object.Z[-1])*0.5]]).ravel()
    z_width = zedges[1:]-zedges[:-1]

    #Average number of DLAs per pixel
    N = z_width*dNdz(object.Z)

    #For a given z, probability of having the density higher than the threshold
    p_nu_z = 1.0-norm.cdf(nu_arr)
    mu = N/p_nu_z

    #Should the "len(skewers)" be the number of skewers or the number of cells in each skewer here?
    #Think it's the number of skewers but will check
    pois = np.random.poisson(mu,size=(object.N_qso,len(mu)))
    dlas = pois*flagged_pixels

    ndlas = np.sum(dlas)
    zdla = np.zeros(ndlas)
    kskw = np.zeros(ndlas)
    dz_dla = np.zeros(ndlas)
    idx = 0

    #In each skewer, look at each potential DLA.
    for nskw,dla in enumerate(dlas):

        #Asess which potential DLA position will be allocated a DLA.
        ind = np.where(dla>0)[0]

        #For each dla, assign it a redshift, a velocity and a column density.
        for ii in ind:
            zdla[idx:idx+dla[ii]] = np.random.uniform(low=(zedges[ii]),high=(zedges[ii+1]),size=dla[ii])
            kskw[idx:idx+dla[ii]] = nskw
            dz_dla[idx:idx+dla[ii]] = object.VEL_rows[nskw,ii]
            idx = idx+dla[ii]

    Ndla = get_N(zdla)
    kskw = kskw.astype('int32')
    MOCKIDs = object.MOCKID[kskw]

    #Make the data into a table HDU
    taux = astropy.table.Table([MOCKIDs,zdla,dz_dla,Ndla],names=('MOCKID','Z_DLA','DZ_DLA','N_HI_DLA'))

    #Only include DLAs where the DLA is at lower z than the QSO
    DLA_Z_QSOs = object.Z_QSO[kskw]
    taux = taux[taux['Z_DLA']<DLA_Z_QSOs]

    object.DLA_table = taux

    return
