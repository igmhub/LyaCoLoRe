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

def get_N(z, Nmin=19.5, Nmax=22, nsamp=100):
    """ Get the column density for a given z"""
    nn = np.linspace(Nmin,Nmax,nsamp)
    probs = dnHD_dz_cumlgN(z,nn).T
    N = np.zeros(len(probs))
    for i in range(0,len(probs)):
        N[i] = np.random.choice(nn,size=1,p=probs[i]/np.sum(probs[i]))
    return N
