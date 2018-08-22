import numpy as np
from astropy.io import fits
import mcfit
import matplotlib.pyplot as plt
import sys

from . import plot_functions

def visual_fit(filename,b_values,beta_values,model,data_parameters,z,compute_b_beta=False):

    mubin_boundaries = [0.0,1.0]
    mubins = []
    for i in range(len(mubin_boundaries)-1):
        mubins += [(mubin_boundaries[i],mubin_boundaries[i+1])]
    #mubins = [(0.0,0.33),(0.33,0.67),(0.67,1.0)]
    mubins = [(0.0,0.33),(0.33,0.67),(0.67,1.0)]
    #find a more sophisticated way to do this
    colours = ['r',(0.5,0.5,0.5),'b']

    #mubins = [(0.0,1.0)]
    #colours = ['r']

    N_bins = len(mubins)

    plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')

    for i,bin in enumerate(mubins):
        mumin = bin[0]
        mumax = bin[1]

        #Wedgize the data
        r,xi_wed,err_wed,cut,_ = plot_functions.get_plot_data(mumin,mumax,filename)

        data_label = 'data, {}<$\mu$<{}'.format(mumin,mumax)
        plt.errorbar(r[cut],xi_wed[cut]*(r[cut]**2),yerr=err_wed[cut]*(r[cut]**2),fmt='o',label=data_label,c=colours[i])

    quantity1 = data_parameters['quantities'][0]
    quantity2 = data_parameters['quantities'][1]

    #TODO: implement compute_b_beta. Need to have different flags for each quantity?
    for b1 in b_values[quantity1]:
        for beta1 in beta_values[quantity1]:
            for b2 in b_values[quantity2]:
                for beta2 in beta_values[quantity2]:

                    r_model,xi_model_values = get_model_xi(model,[b1,b2],[beta1,beta2],data_parameters,z,mubins,b_from_z=False)

                    for i,key in enumerate(xi_model_values.keys()):
                        #model_label = 'b_{}={}, beta_{}={}, b_{}={}, beta_{}={}, mu={}'.format(quantity1,b1,quantity1,beta1,quantity2,b2,quantity2,beta2,key)
                        model_label = 'model, $\mu$={}'.format(key)
                        plt.plot(r_model,xi_model_values[key]*(r_model**2),label=model_label,c=colours[i])
                        #plt.plot(r[cut],(xi_wed[cut])/(np.interp(r[cut],r_model,xi_model_values[key])),label='(RATIO): '+model_label)

    plt.axhline(y=0,color='gray',ls=':')
    plt.xlabel(r'$r / Mpc h^{-1}$')
    plt.ylabel(r'$r^2 \xi(r)$')
    plt.grid()
    plt.legend()
    plt.xlim(0,200)
    #Make this more sensible filename
    plt.savefig('fit_picca_visual.pdf')
    plt.show()

    return

def get_model_xi(model,bs,betas,data_parameters,z,mubins,b_from_z=False):

    b1 = bs[0]
    b2 = bs[1]
    beta1 = betas[0]
    beta2 = betas[1]

    quantity1 = data_parameters['quantities'][0]
    quantity2 = data_parameters['quantities'][1]

    if model == 'no_beta':
        sr = data_parameters['sr']

        #Open xi data files
        file_location = '/global/homes/j/jfarr/Projects/LyaCoLoRe/camb_scripts/camb_xi_10.txt'
        r = np.loadtxt(file_location)[:,0]
        xi = np.loadtxt(file_location)[:,1]

        #Calculate the appropriate scaling
        scaling = get_growth_factor_scaling(z,quantity1)*get_growth_factor_scaling(z,quantity2)
        if b_from_z:
            scaling *= get_bias(z,quantity1)*get_bias(z,quantity2)
        else:
            scaling *= b1*b2

        xi *= scaling

    elif model == 'Slosar11':

        Pk_location = '/global/homes/j/jfarr/Projects/run_CoLoRe/input_files/Pk_CAMB_test.dat'

        Pk_CAMB = np.loadtxt(Pk_location)

        k_old = Pk_CAMB[:,0]
        P_old = Pk_CAMB[:,1]

        #Pick the mu values from the centre of each bin
        mu_values = []
        for mubin in mubins:
            mu_values += [(mubin[0]+mubin[1])/2.]

        k_min = -4
        k_max = 3
        k_num = 5

        xi_values = {}

        for mu in mu_values:

            P_mu_0 = np.polynomial.legendre.legval(mu,[1])
            P_mu_2 = np.polynomial.legendre.legval(mu,[0,0,1])
            P_mu_4 = np.polynomial.legendre.legval(mu,[0,0,0,0,1])

            filename = '/global/homes/j/jfarr/Projects/PhD/xil/xil_{}_{}_{}.txt'.format(k_min,k_max,k_num)
            data = np.loadtxt(filename)

            r = data[:,0]
            xi0 = data[:,1]
            xi2 = data[:,2]
            xi4 = data[:,3]

            C0 = get_C0(beta1,beta2)
            C2 = get_C2(beta1,beta2)
            C4 = get_C4(beta1,beta2)

            scaling = get_growth_factor_scaling(z,quantity1)*get_growth_factor_scaling(z,quantity2)

            if b_from_z:
                scaling *= get_bias(z,quantity1)*get_bias(z,quantity2)
            else:
                scaling *= b1*b2

            xi = scaling*(C0*xi0*P_mu_0 + C2*xi2*P_mu_2 + C4*xi4*P_mu_4)

            new_xi_value = {mu: xi}
            xi_values = {**xi_values,**new_xi_value}

    return r, xi_values

def get_C0(B1,B2):
    return 1 + (1/3)*(B1+B2) + (1/5)*B1*B2

def get_C2(B1,B2):
    return (2/3)*(B1+B2) + (4/7)*B1*B2

def get_C4(B1,B2):
    return (8/35)*B1*B2

def get_growth_factor_scaling(z,quantity,location=None):

    if location == None:
        location = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/test/'
    if quantity == 'G':
        D_at_z0 = 1
        D_at_zval = 1
    elif quantity in ['D','F','q']:
        h = fits.open('/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/test/master.fits')
        D = h[2].data['D']
        z_D = h[2].data['Z']
        D_at_z0 = np.interp(0,z_D,D)
        D_at_zval = np.interp(z,z_D,D)
    else:
        print('quantity not recognised')

    return D_at_zval/D_at_z0

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
