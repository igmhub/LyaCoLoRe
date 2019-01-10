from __future__ import print_function, division
import numpy as np
#import astropy.io.fits as fits
from scipy.stats import norm
from scipy.interpolate import interp1d, interp2d
import astropy.table
import os
import fitsio
import glob
#import matplotlib.pyplot as plt
import argparse

try:
    from pyigm.fN.fnmodel import FNModel
    fN_default = FNModel.default_model()
    fN_default.zmnx = (0.,4)
    fN_cosmo = fN_default.cosmo
    use_pyigm = True
except:
    use_pyigm = False

def doloop(dlas_in_cell,velocity,zedges, dla_z, dla_skw_id, dla_rsd_dz, dla_count):
    """ Auxiliary function to perform the loop to populate the DLA cells"""
    for skw_id,dla in enumerate(dlas_in_cell):
        #Find cells that will be allocated at least one DLA
        dla_cells = np.where(dla>0)[0]
        #For each dla, assign it a redshift, a velocity and a column density.
        for cell in dla_cells:
            dla_z[dla_count:dla_count+dla[cell]] = zedges[cell]+(zedges[cell+1]-zedges[cell])*np.random.random(size=dla[cell])
            dla_skw_id[dla_count:dla_count+dla[cell]] = skw_id
            dla_rsd_dz[dla_count:dla_count+dla[cell]] = velocity[skw_id,cell]
            dla_count = dla_count+dla[cell]
    return dla_z, dla_skw_id, dla_rsd_dz, dla_count

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

def flag_DLA(zq,z_cells,deltas,nu_arr,sigma_g,zlow):
    """ Flag the pixels in a skewer where DLAs are possible"""
    # find cells with density above threshold
    flag = deltas > nu_arr*sigma_g
    # mask cells with z > z_qso, where DLAs would not be observed
    Nq=len(zq)
    for i in range(Nq):
        low_z = (z_cells < zq[i]) & (z_cells > zlow)
        flag[i,:] *= low_z
    return flag

#number per unit redshift from minimum lg(N) in file (17.2) to argument
# Reading file from https://arxiv.org/pdf/astro-ph/0407378.pdf

def dnHD_dz_cumlgN(z,logN):
    tab = astropy.table.Table.read(os.path.abspath('LyaCoLoRe/example_data/zheng_cumulative.overz'),format='ascii')
    y = interp2d(tab['col1'],tab['col2'],tab['col3'],fill_value=None)
    return y(z,logN)

def dNdz(z, Nmin=20.0, Nmax=22.5):
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


def get_N(z, Nmin=20.0, Nmax=22.5, nsamp=100):
    """ Get random column densities for a given z
    """

    # number of DLAs we want to generate
    Nz = len(z)
    nn = np.linspace(Nmin,Nmax,nsamp)
    probs = np.zeros([Nz,nsamp])
    if use_pyigm:
        auxfN = fN_default.evaluate(nn,z)
        #auxfN = (np.cumsum(10**auxfN, axis=0)/np.sum(10**auxfN, axis=0)).T
        probs = (np.exp(auxfN)/np.sum(np.exp(auxfN), axis=0)).T
        #plt.plot(nn,auxfN.T)
    else:
        probs_low = dnHD_dz_cumlgN(z,nn[:-1]).T
        probs_high = dnHD_dz_cumlgN(z,nn[1:]).T
        probs[:,1:] = probs_high-probs_low
    NHI = np.zeros(Nz)
    for i in range(Nz):
        #if use_pyigm:
        #    nfunc = interp1d(auxfN[i],nn,fill_value='extrapolate')
        #    NHI[i] = nfunc(np.random.uniform())
        #else:
        #    NHI[i] = np.random.choice(nn,size=1,p=probs[i]/np.sum(probs[i]))+(nn[1]-nn[0])*np.random.random(size=1)
        NHI[i] = np.random.choice(nn,size=1,p=probs[i]/np.sum(probs[i]))+(nn[1]-nn[0])*np.random.random(size=1)
    return NHI


def add_DLA_table_to_object_Saclay(fname,fname_cosmo,fname_sigma,dNdz_arr,dla_bias=20.0,extrapolate_z_down=None,Nmin=20.0,Nmax=22.5,seed=123,zlow=1.8):
    np.random.seed(seed)
    hdulist = fitsio.FITS(fname) # Open the file
    qso = hdulist[1].read() # Read the QSO table
    lam = hdulist[2].read() # Read the vector with the wavelenghts corresponding to each cell
    #velocity = hdulist[3].read() # Placeholder
    cosmo_hdu = fitsio.FITS(fname_cosmo)[1].read_header() # Reading the cosmological parameters used for the simulation
    #Quasar redshift for each skewer
    zq = qso['Z']
    #Redshift of each cell in the skewer
    z_cell = lam / 1215.67 - 1
    Oc = cosmo_hdu['OM']-cosmo_hdu['OB'] # Omega_c
    Ob = cosmo_hdu['OB'] # Omega_b
    h = cosmo_hdu['H'] # h
    Ok = cosmo_hdu['OK'] # Omega_k
    #Linear growth rate of each cell in the skewer
    D_cell = hdulist['GROWTHF'].read()   
    #Setup bias as a function of redshift
    y = interp1d(z_cell,D_cell)
    bias = dla_bias/(D_cell)*y(2.25)
    sigma_g = fitsio.FITS(fname_sigma)[0].read_header()['SIGMA']
    nu_arr = nu_of_bD(bias*D_cell)
    #Figure out cells that could host a DLA, based on Gaussian fluctuation
    deltas = hdulist[4].read()
    velocity = np.zeros_like(deltas)
    hdulist.close()
    flagged_cells = flag_DLA(zq,z_cell,deltas,nu_arr,sigma_g,zlow)
    flagged_cells[deltas==-1e6]=False
    #Edges of the z bins
    if extrapolate_z_down and extrapolate_z_down<z_cell[0]:
        zedges = np.concatenate([[extrapolate_z_down],(z_cell[1:]+z_cell[:-1])*0.5,[z_cell[-1]+(-z_cell[-2]+z_cell[-1])*0.5]]).ravel()
    else:
        zedges = np.concatenate([[z_cell[0]],(z_cell[1:]+z_cell[:-1])*0.5,[z_cell[-1]+(-z_cell[-2]+z_cell[-1])*0.5]]).ravel()
    z_width = zedges[1:]-zedges[:-1]

    #Get the average number of DLAs per cell, from the column density dist.
    mean_N_per_cell = z_width*dNdz_arr
    #For a given z, probability of having the density higher than the threshold
    p_nu_z = 1.0-norm.cdf(nu_arr)
    #Define mean of the Poisson distribution (per cell)
    #mu = mean_N_per_cell/p_nu_z
    mu = mean_N_per_cell*(1+bias*deltas)
    mu[~flagged_cells]=0
    #Select cells that will hold a DLA, drawing from the Poisson distribution
    pois = np.random.poisson(mu)#,size=(len(zq),len(mu)))
    #Number of DLAs in each cell (mostly 0, several 1, not many with >1)
    dlas_in_cell = pois*flagged_cells
    ndlas = np.sum(dlas_in_cell)
    #Store information for each of the DLAs that will be added
    dla_z = np.zeros(ndlas)
    dla_skw_id = np.zeros(ndlas, dtype='int32')
    dla_rsd_dz = np.zeros(ndlas)
    dla_count = 0
    dla_z, dla_skw_id, dla_rsd_dz, dla_count = doloop(dlas_in_cell, velocity, zedges, dla_z, dla_skw_id, dla_rsd_dz, dla_count)
    dla_NHI = get_N(dla_z,Nmin=Nmin,Nmax=Nmax)

    #global id for the skewers
    MOCKIDs = qso['THING_ID'][dla_skw_id]
    ZQSO = zq[dla_skw_id]
    #Make the data into a table HDU
    dla_table = astropy.table.Table([MOCKIDs,dla_z,dla_rsd_dz,dla_NHI,ZQSO],names=('MOCKID','Z_DLA','DZ_DLA','N_HI_DLA','Z_QSO'))

    return dla_table

######

# Options and main

######

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input_path', type = str, default = None, required = True,
                    help='Path to input directory tree to explore, e.g., /global/cscratch1/sd/*/spectra/*')
parser.add_argument('--output_file', type = str, default = None, required = True,
                    help='Output filename')
parser.add_argument('--input_pattern', type = str, default = 'spectra_merged*.fits',
                    help='Filename pattern')
parser.add_argument('--fname_cosmo', type = str, default = None, required = True,
                    help='Path to file with cosmological parameter information')
parser.add_argument('--fname_sigma', type = str, default = None, required = True,
                    help='Path to file with information sigma(0) (initial density field RMS)')
parser.add_argument('--nmin', type = float, default=20,
                    help='Minimum value of log(NHI) to consider')
parser.add_argument('--nmax', type = float, default=22.5,
                    help='Maximum value of log(NHI) to consider')
parser.add_argument('--dla_bias', type = float, default=2.0,
                    help='DLA bias at z=2.25')
args = parser.parse_args()

flist = glob.glob(os.path.join(args.input_path,args.input_pattern))
print('Will read', len(flist),' files')
hdulist = fitsio.FITS(flist[0])
lam = hdulist[2].read()
cosmo_hdu = fitsio.FITS(args.fname_cosmo)[1].read_header()
z_cell = lam / 1215.67 - 1.
dNdz_arr = dNdz(z_cell, Nmin=args.nmin, Nmax=args.nmax)

for i, fname in enumerate(flist):
    aux = add_DLA_table_to_object_Saclay(fname, args.fname_cosmo, args.fname_sigma, dNdz_arr, args.dla_bias) 
    if i==0:
        out_table = aux
    else:
        out_table = astropy.table.vstack([out_table, aux])
    if i%500==0:
        print('Read %d of %d' %(i,len(flist)))

out_table.write(args.output_file, overwrite=True)
