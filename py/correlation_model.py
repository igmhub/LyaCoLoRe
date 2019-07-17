import numpy as np
from astropy.io import fits
import mcfit
import matplotlib.pyplot as plt
import sys

def get_model_xi(model,q1,q2,bias1,bias2,beta1,beta2,z,mubin,sr=0.0):

    if model == 'no_beta':

        #Open xi data files
        file_location = '/global/homes/j/jfarr/Projects/LyaCoLoRe/camb_scripts/camb_xi_10.txt'
        r = np.loadtxt(file_location)[:,0]
        xi = np.loadtxt(file_location)[:,1]

        #Calculate the appropriate scaling
        scaling = get_growth_factor_scaling(z,q1)*get_growth_factor_scaling(z,q2)
        if b_from_z:
            scaling *= get_bias(z,q1)*get_bias(z,q2)
        else:
            scaling *= b1*b2

        xi *= scaling

    elif model == 'Slosar11':

        #Pick the mu values from the centre of each bin
        mu = (mubin[0]+mubin[1])/2.

        k_min = -4
        k_max = 3
        k_num = 5

        P_mu_0 = np.polynomial.legendre.legval(mu,[1])
        P_mu_2 = np.polynomial.legendre.legval(mu,[0,0,1])
        P_mu_4 = np.polynomial.legendre.legval(mu,[0,0,0,0,1])

        # TODO: bring these closer: atm it's useless for anyone other than me
        filename = '/global/homes/j/jfarr/Projects/PhD/xil/xil_{}_{}_{}.txt'.format(k_min,k_max,k_num)
        data = np.loadtxt(filename)

        r = data[:,0]
        xi0 = data[:,1]
        xi2 = data[:,2]
        xi4 = data[:,3]

        C0 = get_C0(beta1,beta2)
        C2 = get_C2(beta1,beta2)
        C4 = get_C4(beta1,beta2)

        scaling = get_growth_factor_scaling(z,q1)*get_growth_factor_scaling(z,q2)
        print(get_growth_factor_scaling(z,q1))
        scaling *= bias1*bias2

        xi = scaling*(C0*xi0*P_mu_0 + C2*xi2*P_mu_2 + C4*xi4*P_mu_4)

    return r, xi

def get_model_xi_grid(model,q1,q2,bias1,bias2,beta1,beta2,z,rp_grid,rt_grid,sr=0.0):

    if model == 'no_beta':

        #Open xi data files
        file_location = './camb_scripts/camb_xi_10.txt'
        r = np.loadtxt(file_location)[:,0]
        xi = np.loadtxt(file_location)[:,1]

        #Calculate the appropriate scaling
        scaling = get_growth_factor_scaling(z,q1)*get_growth_factor_scaling(z,q2)
        if b_from_z:
            scaling *= get_bias(z,q1)*get_bias(z,q2)
        else:
            scaling *= b1*b2

        xi *= scaling

    elif model == 'Slosar11':

        r_grid = np.sqrt(rp_grid**2 + rt_grid**2)
        mu_grid = rp_grid / r_grid

        k_min = -4
        k_max = 3
        k_num = 5

        P_mu_0 = np.polynomial.legendre.legval(mu_grid,[1])
        P_mu_2 = np.polynomial.legendre.legval(mu_grid,[0,0,1])
        P_mu_4 = np.polynomial.legendre.legval(mu_grid,[0,0,0,0,1])

        # TODO: bring these closer: atm it's useless for anyone other than me
        filename = '../PhD/xil/xil_{}_{}_{}.txt'.format(k_min,k_max,k_num)
        data = np.loadtxt(filename)

        r = data[:,0]
        xi0 = data[:,1]
        xi2 = data[:,2]
        xi4 = data[:,3]

        xi0_grid = np.interp(r_grid,r,xi0)
        xi2_grid = np.interp(r_grid,r,xi2)
        xi4_grid = np.interp(r_grid,r,xi4)

        C0 = get_C0(beta1,beta2)
        C2 = get_C2(beta1,beta2)
        C4 = get_C4(beta1,beta2)

        scaling = get_growth_factor_scaling(z,q1)*get_growth_factor_scaling(z,q2)
        scaling *= bias1*bias2

        xi_grid = scaling*(C0*xi0_grid*P_mu_0 + C2*xi2_grid*P_mu_2 + C4*xi4_grid*P_mu_4)

    return r_grid, xi_grid



def get_C0(B1,B2):
    return 1 + (1/3)*(B1+B2) + (1/5)*B1*B2

def get_C2(B1,B2):
    return (2/3)*(B1+B2) + (4/7)*B1*B2

def get_C4(B1,B2):
    return (8/35)*B1*B2

def get_growth_factor_scaling(z,quantity,location=None):

    if location == None:
        location = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v6/v6.0.0/'
    if quantity in ['G']:
        D_at_z0 = 1
        D_at_zval = 1
    elif quantity in ['D','F','T','q']:
        h = fits.open(location+'/master.fits')
        D = h['COSMO_COL'].data['D']
        z_D = h['COSMO_COL'].data['Z']
        D_at_z0 = np.interp(0,z_D,D)
        D_at_zval = np.interp(z,z_D,D)
    else:
        print('quantity not recognised')

    return (D_at_zval/D_at_z0)






def get_bias(z,quantity):

    if quantity == 'q':
        bias_data = np.loadtxt('/global/homes/j/jfarr/Projects/run_CoLoRe/input_files/Bz_qso_G18.txt')
        z_bq = bias_data[:,0]
        bq = bias_data[:,1]
        bq_at_zval = np.interp(z,z_bq,bq)
    elif quantity in ['G','D']:
        bq_at_zval = 1
    else:
        print('quantity not recognised')

    return bq_at_zval

def beta_QSO_kaiser(z,b,Om_z0=0.3147):

    Om = Om_z0 * ((1+z)**3) / (Om_z0 * ((1+z)**3) + 1 - Om_z0)
    Ol = 1 - Om
    # TODO: replace this with the analytic expression using D from file
    f = (Om**0.6) + (Ol/70.)*(1 + Om/2.) #using review eqn 4.14
    beta = f/b

    return beta
