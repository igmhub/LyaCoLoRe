from __future__ import print_function, division
import numpy as np
import astropy.io.fits as fits
from scipy.stats import norm
from scipy.interpolate import interp1d, interp2d
from optparse import OptionParser
import astropy.table
import os
import matplotlib.pyplot as plt

parser  = OptionParser()

parser.add_option("--dla_bias", dest="dla_bias", type=float, default=2.0,
    help = "Constant multiplied to the inverse growth to get the bias")
parser.add_option("--input-file", dest="input_file", type=str, default=None,
    help = "Input path")
parser.add_option("--output-file", dest="output_name", type=str, default=None,
    help = "Output path, if not passed it will overwrite the input file")

(o, args) = parser.parse_args()

def nu_of_bD(b):
    """ Compute the Gaussian field threshold for a given bias"""
    nu = np.linspace(-10,100,500) # Generous range to interpolate
    p_nu = norm.pdf(nu)
    galaxy_mean = 1.0-norm.cdf(nu)
    b_nu = np.zeros(nu.shape)
    indices = galaxy_mean!=0
    b_nu[indices] = p_nu[indices]/galaxy_mean[indices]
    #plt.plot(b_nu[indices],nu[indices]);plt.grid();plt.show()

    y = interp1d(b_nu[indices],nu[indices])
    #y = np.interp(b,b_nu[indices],nu[indices])
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

if o.output_name is not None:
    fileout = o.output_name
else:
    fileout = o.input_file

skewers = fits.open(o.input_file)[2].data
v_skw = fits.open(o.input_file)[3].data
z, bias, D = get_bias_z(o.input_file,o.dla_bias)
sigma_g = get_sigma_g(o.input_file)
nu_arr = nu_of_bD(bias*D)
flagged_pixels = flag_DLA(skewers,nu_arr,sigma_g)
zedges = np.concatenate([[0],(z[1:]+z[:-1])*0.5,[z[-1]+(-z[-2]+z[-1])*0.5]]).ravel()
z_width = zedges[1:]-zedges[:-1]
N = z_width*dNdz(z) # Average number of DLAs per pixel
p_nu_z = 1.0-norm.cdf(nu_arr) # For a given z, probability of having the density higher than the threshold
mu = N/p_nu_z
pois = np.random.poisson(mu,size=(len(skewers),len(mu)))
dlas = pois*flagged_pixels
ndlas = np.sum(dlas)
zdla = np.zeros(ndlas)
kskw = np.zeros(ndlas)
dz_dla = np.zeros(ndlas)
idx = 0
for nskw,dla in enumerate(dlas):
    ind = np.where(dla>0)[0]
    for ii in ind:
        zdla[idx:idx+dla[ii]]=np.random.uniform(low=(zedges[ii]),high=(zedges[ii+1]),size=dla[ii])
        kskw[idx:idx+dla[ii]]=nskw
        dz_dla[idx:idx+dla[ii]]=v_skw[nskw,ii]
        idx = idx+dla[ii]
Ndla = get_N(zdla)
taux = astropy.table.Table([kskw,zdla,dz_dla,Ndla],names=('SKEWER_NUMBER','Z_DLA','DZ_DLA','N_HI_DLA'))
new_hdu = fits.hdu.BinTableHDU(data=taux,name='DLA')
hdulist = fits.open(o.input_file)
hdulist.append(new_hdu)
hdulist.writeto(fileout, overwrite=True)
