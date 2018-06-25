import numpy as np
from astropy.io import fits

import convert

lya = 1215.67

#Function to return the P1D from Palanque-Delabrouille et al. (2013)
#copied from lyaforecast
def P1D_z_kms_PD2013(z,k_kms):
    """Fitting formula for 1D P(z,k) from Palanque-Delabrouille et al. (2013).
        Wavenumbers and power in units of km/s. Corrected to be flat at low-k"""
    # numbers from Palanque-Delabrouille (2013)
    A_F = 0.064
    n_F = -2.55
    alpha_F = -0.1
    B_F = 3.55
    beta_F = -0.28
    k0 = 0.009
    z0 = 3.0
    n_F_z = n_F + beta_F * np.log((1+z)/(1+z0))
    # this function would go to 0 at low k, instead of flat power
    k_min=k0*np.exp((-0.5*n_F_z-1)/alpha_F)
    k_kms = np.fmax(k_kms,k_min)
    exp1 = 3 + n_F_z + alpha_F * np.log(k_kms/k0)
    toret = np.pi * A_F / k0 * pow(k_kms/k0, exp1-1) * pow((1+z)/(1+z0), B_F)
    return toret

#Function to return the mean value of F at a given redshift.
#Equation from F-R2012, equation 2.11
def get_mean_F_model(z):
    mean_F = np.exp((np.log(0.8))*(((1+z)/3.25)**3.2))
    return mean_F

#Function to find the value of alpha required to match mean_F to a specified value.
def find_alpha(sigma_G,mean_F_required,beta,D,alpha_log_low=-3.0,alpha_log_high=10.0,tolerance=0.0001,max_iter=30):
    #print('---> mean_F required={:2.2f}'.format(mean_F_required))
    count = 0
    exit = 0
    while exit == 0 and count < max_iter:
        alpha_log_midpoint = (alpha_log_low + alpha_log_high)/2.0

        mean_F_al,sigma_dF_al = get_flux_stats(sigma_G,10**alpha_log_low,beta,D,mean_only=True)
        mean_F_am,sigma_dF_am = get_flux_stats(sigma_G,10**alpha_log_midpoint,beta,D,mean_only=True)
        mean_F_ah,sigma_dF_ah = get_flux_stats(sigma_G,10**alpha_log_high,beta,D,mean_only=True)

        #print('---> alphas=({:2.2f},{:2.2f},{:2.2f}) gives mean_F=({:2.2f},{:2.2f},{:2.2f})'.format(10**alpha_log_low,10**alpha_log_midpoint,10**alpha_log_high,mean_F_al,mean_F_am,mean_F_ah))

        if np.sign(mean_F_al-mean_F_required) * np.sign(mean_F_am-mean_F_required) > 0:
            alpha_log_low = alpha_log_midpoint
        else:
            alpha_log_high = alpha_log_midpoint

        if abs(mean_F_am/mean_F_required - 1) < tolerance:
            exit = 1
        else:
            count += 1

    if exit == 0:
        # TODO: something other than print here. Maybe make a log of some kind?
        print('\nvalue of mean_F did not converge to within tolerance: error is {:3.2%}'.format(mean_F_am/mean_F_required - 1))

    alpha = 10**alpha_log_midpoint
    mean_F,sigma_dF = get_flux_stats(sigma_G,alpha,beta,D)

    return alpha,mean_F,sigma_dF

#Function to find the values of alpha and sigma_G required to match mean_F and sigma_dF to specified values.
def find_sigma_G(mean_F_required,sigma_dF_required,beta,D,sigma_G_start=0.001,step_size=1.0,tolerance=0.001,max_steps=30):
    #print('sigma_dF required={:2.2f}'.format(sigma_dF_required))
    #print(' ')
    #print(' ')
    count = 0
    exit = 0
    sigma_G = sigma_G_start

    while exit == 0 and count < max_steps:

        alpha,mean_F,sigma_dF = find_alpha(sigma_G,mean_F_required,beta,D)

        if abs(sigma_dF/sigma_dF_required - 1) < tolerance:
            exit = 1
            #print('sigma_G={:2.4f} gives sigma_dF={:2.4f}. Satisfied. Exiting...'.format(sigma_G,sigma_dF))
        elif sigma_dF < sigma_dF_required:
            #print('sigma_G={:2.4f} gives sigma_dF={:2.4f}. Too low. Stepping forwards...'.format(sigma_G,sigma_dF))
            sigma_G += step_size
            count += 1
        elif sigma_dF > sigma_dF_required:
            #print('sigma_G={:2.4f} gives sigma_dF={:2.4f}. Too high. Stepping backwards...'.format(sigma_G,sigma_dF))
            sigma_G_too_high = sigma_G
            sigma_G -= step_size
            step_size = step_size/10.0
            sigma_G += step_size
            count += 1

        #print('error: ',(sigma_dF/sigma_dF_required - 1))

        """
        sigma_G_log_low = -3.0
        sigma_G_log_high = 1.0

        sigma_G_log_midpoint = (sigma_G_log_low + sigma_G_log_high)/2.0

        alpha_sGl,mean_F_sGl,sigma_dF_sGl = find_alpha(10**sigma_G_log_low,mean_F_required,beta,D)
        alpha_sGm,mean_F_sGm,sigma_dF_sGm = find_alpha(10**sigma_G_log_midpoint,mean_F_required,beta,D)
        alpha_sGh,mean_F_sGh,sigma_dF_sGh = find_alpha(10**sigma_G_log_high,mean_F_required,beta,D)

        print('sigma_Gs=({:2.2f},{:2.2f},{:2.2f}) gives sigma_dFs=({:2.2f},{:2.2f},{:2.2f})'.format(10**sigma_G_log_low,10**sigma_G_log_midpoint,10**sigma_G_log_high,sigma_dF_sGl,sigma_dF_sGm,sigma_dF_sGh))

        if np.sign(sigma_dF_sGl-sigma_dF_required) * np.sign(sigma_dF_sGm-sigma_dF_required) > 0:
            sigma_G_log_low = sigma_G_log_midpoint
        else:
            sigma_G_log_high = sigma_G_log_midpoint

        if abs(sigma_dF_sGm/sigma_dF_required - 1) < tolerance:
            exit = 1
        else:
            count += 1
        """

    #print('Testing finished. Final values are:')
    #print('sigma_G={:2.4f} gives sigma_dF={:2.4f}.'.format(sigma_G,sigma_dF))
    #print('error: ',(sigma_dF/sigma_dF_required - 1))

    if exit == 0:
        # TODO: something other than print here. Maybe make a log of some kind?
        print('\nvalue of sigma_dF did not converge to within tolerance: error is {:3.2%}'.format(sigma_dF/sigma_dF_required - 1))
        sigma_G = (sigma_G+sigma_G_too_high)/2.0
        alpha,mean_F,sigma_dF = find_alpha(sigma_G,mean_F_required,beta,D)

    #print('Final check finished. Final values are:')
    #print('sigma_G={:2.4f} gives sigma_dF={:2.4f}.'.format(sigma_G,sigma_dF))
    #print('error: ',(sigma_dF/sigma_dF_required - 1))
    #print(' ')
    """
    alpha = alpha_sGm
    sigma_G = 10**sigma_G_log_midpoint
    mean_F = mean_F_sGm
    sigma_dF = sigma_dF_sGm
    """
    return alpha,sigma_G,mean_F,sigma_dF

#Function to calculate mean_F and sigma_dF for given values of sigma_G, alpha and beta.
def get_flux_stats(sigma_G,alpha,beta,D,mean_only=False,int_lim_fac=10.0):

    int_lim = sigma_G*int_lim_fac

    delta_G_integral = np.linspace(-int_lim,int_lim,10**4)
    delta_G_integral = np.reshape(delta_G_integral,(1,delta_G_integral.shape[0]))

    prob_delta_G = (1/((np.sqrt(2*np.pi))*sigma_G))*np.exp(-(delta_G_integral**2)/(2*(sigma_G**2)))

    density_integral = convert.gaussian_to_lognormal_delta(delta_G_integral,sigma_G,D) + 1
    F_integral = convert.density_to_flux(density_integral,alpha,beta)

    mean_F = np.trapz(prob_delta_G*F_integral,delta_G_integral)[0]

    if mean_only == False:
        delta_F_integral = F_integral/mean_F - 1
        integrand = prob_delta_G*(delta_F_integral**2)
        sigma_dF = (np.sqrt(np.trapz(integrand,delta_G_integral)[0]))
    else:
        sigma_dF = None

    return mean_F, sigma_dF

#Function to integrate under the 1D power spectrum to return the value of sigma_dF at a given redshift.
def get_sigma_dF_P1D(z,l_hMpc=0.25,Om=0.3):
    #Choose log spaced values of k
    k_hMpc_max = 100.0/l_hMpc
    k_hMpc = np.logspace(-5,np.log10(k_hMpc_max),10**5)

    # TODO: generalise the conversion in here
    # need to go from Mpc/h to km/s, using dv / dX = H(z) / (1+z)
    # we will define H(z) = 100 h E(z)
    # with E(z) = sqrt(Omega_m(1+z)^3 + Omega_L), and assume flat universe
    E_z = np.sqrt(Om*(1+z)**3 + (1-Om))
    dkms_dhMpc = 100. * E_z / (1+z)

    # transform input wavenumbers to s/km
    k_kms = k_hMpc / dkms_dhMpc

    # get power in units of km/s
    pk_kms = P1D_z_kms_PD2013(z,k_kms)

    # transform to h/Mpc
    pk_hMpc = pk_kms / dkms_dhMpc

    # compute Fourier transform of Top-Hat filter of size l_hMpc
    W_hMpc = np.sinc((k_hMpc*l_hMpc)/(2*np.pi))
    sigma_dF = np.sqrt((1/np.pi)*np.trapz((W_hMpc**2)*pk_hMpc,k_hMpc))

    return sigma_dF
