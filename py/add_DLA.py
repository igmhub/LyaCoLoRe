from __future__ import print_function, division
import numpy as np
import astropy.io.fits as fits
from scipy.stats import norm
from scipy.interpolate import interp1d
from optparse import OptionParser

parser  = OptionParser()

parser.add_option("--dla_bias", dest="dla_bias", type=float, default=2.0,
    help = "Constant multiplied to the inverse growth to get the bias")
parser.add_option("--input-file", dest="input_file", type=str, default=None,
    help = "Input path")
parser.add_option("--output-file", dest="output_file", type=str, default=None,
    help = "Output path, if not passed it will overwrite the input file")

(o, args) = parser.parse_args()

def nu_of_b(b):
    """ Compute the Gaussian field threshold for a given bias"""
    nu = np.linspace(-10,100,500) # Generous range to interpolate
    p_nu = norm.pdf(nu)
    galaxy_mean = 1.0-norm.cdf(nu)
    b_nu = p_nu / galaxy_mean
    y = interp1d(b_nu,nu)
    return y(b)

def get_bias_z(fname,dla_bias):
    """ Given a path, read the z array there and return a bias inversely
    proportional to the growth"""
    cosmo = fits.open(fname)[4].data
    z = cosmo['Z']
    D = cosmo['D']
    bias = dla_bias/D
    return z, bias

def flag_DLA(skewer,nu_arr):
    """ Flag the pixels in a skewer where DLAs are possible"""
    flag = skewer > nu_arr
    return flag

if o.output_name is not None:
    fileout = o.output_name
else:
    fileout = o.input_name

skewers = fits.open(o.input_file)[2].data
z, bias = get_bias_z(o.input_file,o.dla_bias)
nu_arr = nu_of_b(bias)
flagged_pixels = flag_DLA(skewer,nu_arr)
new_hdu = fits.hdu.ImageHDU(flagged_pixels.astype(np.uint8),name='DLA_FLAGS')
hdulist = fits.open(o.input_name)
hdulist.append(new_hdu)
hdulist.writeto(fileout)
