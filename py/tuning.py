import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

import convert
import Pk1D
import general

lya = 1215.67

################################################################################
"""
Below: new tuning, measurement based
"""


class tuning_parameters:

    def __init__(self,):
        return

    @classmethod
    def add_parameter(name,value):



        return cls()


class measurement:
    def __init__(self,parameter_ID,z_value,z_width,N_skewers,n,k1,alpha,beta,sigma_G,pixels=[],mean_F=None,k_kms=None,Pk_kms=None,cf=None):
        self.parameter_ID = parameter_ID
        self.z_value = z_value
        self.z_width = z_width
        self.N_skewers = N_skewers

        self.n = n
        self.k1 = k1
        self.alpha = alpha
        self.beta = beta
        self.sigma_G = sigma_G

        self.pixels = pixels

        self.mean_F = mean_F
        self.k_kms = k_kms
        self.Pk_kms = Pk_kms
        self.cf = cf
        return
    def add_Pk1D_measurement(self,pixel_object):
        F_rows = pixel_object.F_rows
        mean_F = np.average(F_rows)
        delta_F_rows = F_rows/mean_F - 1
        IVAR_rows = pixel_object.IVAR_rows
        R_hMpc = pixel_object.R
        z = pixel_object.Z
        k_kms, Pk_kms, var_kms = Pk1D.get_Pk1D(delta_F_rows,IVAR_rows,R_hMpc,z,self.z_value,z_width=self.z_width)
        self.k_kms = k_kms
        self.Pk_kms = Pk_kms
        return
    def add_mean_F_measurement(self,pixel_object):
        self.mean_F = pixel_object.get_mean_flux(z_value=self.z_value,z_width=self.z_width)
        return
    def add_Pk1D_chi2(self,min_k=None,max_k=None,denom="krange10"):
        model_Pk_kms = P1D_z_kms_PD2013(self.k_kms,self.z_value)
        if min_k:
            min_j = max(np.searchsorted(self.k_kms,min_k) - 1,0)
        else:
            min_j = 0
        if max_k:
            max_j = np.searchsorted(self.k_kms,max_k)
        else:
            max_j = -1
        if denom == "uniform5":
            denom = (0.05 * model_Pk_kms)**2
        elif denom == "uniform10":
            denom = (0.10 * model_Pk_kms)**2
        elif denom == "krange5":
            eps = 10**6 * np.ones(self.k_kms.shape)
            eps[min_j:max_j] *= 0.05 / 10**6
            denom = (eps * model_Pk_kms)**2
        elif denom == "krange10":
            eps = 10**6 * np.ones(self.k_kms.shape)
            eps[min_j:max_j] *= 0.1 / 10**6
            denom = (eps * model_Pk_kms)**2
        chi2 = np.sum(((self.Pk_kms - model_Pk_kms)**2)/denom)
        self.Pk_kms_chi2 = chi2
        return
    def add_mean_F_chi2(self,min_k=None,max_k=None,eps=0.1):
        model_mean_F = get_mean_F_model(self.z_value)
        denom = (eps * model_mean_F)**2
        chi2 = np.sum(((self.mean_F - model_mean_F)**2)/denom)
        self.mean_F_chi2 = chi2
        return
    def add_total_chi2(self):
        chi2 = self.Pk_kms_chi2 + self.mean_F_chi2
        self.total_chi2 = chi2
        return
    @classmethod
    def combine_measurements(cls,m1,m2):
        if general.confirm_identical(m1.parameter_ID,m2.parameter_ID,item_name='parameter_ID'):
            parameter_ID = m1.parameter_ID
            n = m1.n
            k1 = m1.k1
            alpha = m1.alpha
            beta = m1.beta
            sigma_G = m1.sigma_G
        if general.confirm_identical(m1.z_value,m2.z_value,item_name='z_value'):
            z_value = m1.z_value
        if general.confirm_identical(m1.z_width,m2.z_width,item_name='z_width'):
            z_width = m1.z_width
        N_skewers = m1.N_skewers + m2.N_skewers
        pixels = m1.pixels + m2.pixels
        mean_F = (m1.mean_F*m1.N_skewers + m2.mean_F*m2.N_skewers)/(m1.N_skewers + m2.N_skewers)
        k_kms_check = True
        if general.confirm_identical(m1.k_kms,m2.k_kms,item_name='k_kms',array=True):
            k_kms = m1.k_kms
        Pk_kms = (m1.Pk_kms*m1.N_skewers + m2.Pk_kms*m2.N_skewers)/(m1.N_skewers + m2.N_skewers)
        #May need to work on this?
        cf = None
        return measurement(parameter_ID,z_value,z_width,N_skewers,n,k1,alpha,beta,sigma_G,pixels=pixels,mean_F=mean_F,k_kms=k_kms,Pk_kms=Pk_kms,cf=cf)
    """
    def add_Pk1D_error(self,max_k=None):
        model_Pk_kms = P1D_z_kms_PD2013(self.k_kms,self.z_value)
        if max_k:
            max_j = np.searchsorted(self.k_kms,max_k)
        else:
            max_j = -1
        Pk_kms_errors = (self.Pk_kms[:max_j] - model_Pk_kms[:max_j])/model_Pk_kms[:max_j]
        average_error = np.average(Pk_kms_errors)
        self.Pk_kms_error = average_error
        return
    def add_mean_F_error(self):
        model_mean_F = get_mean_F_model(self.z_value)
        mean_F_error = (self.mean_F - model_mean_F)/model_mean_F
        self.mean_F_error = mean_F_error
        return
    """

class measurement_set:
    def __init__(self,measurements=[]):
        self.measurements = measurements
        """
        self.pixels = []
        self.IDs = []
        self.z_values = []
        measurement_index = []
        for m in self.measurements:
            self.pixels += m.pixels
            self.IDs
            measurement_index += [(m.z_value,m.pixels,m.parameter_ID)]
        self.measurement_index = np.array(measurement_index,dtype=dtype)
        """
        return
    def add_measurement(self,measurement):
        self.measurements += [measurement]
        return
    def z_filter(self,z_value):
        filtered_measurements = []
        for m in self.measurements:
            if m.z_value == z_value:
                filtered_measurements += [m]
        return measurement_set(filtered_measurements)
    def s_filter(self,n,k1):
        filtered_measurements = []
        for m in self.measurements:
            if m.n == n and m.k1 == k1:
                filtered_measurements += [m]
        return measurement_set(filtered_measurements)
    def t_filter(self,alpha,beta,sigma_G):
        filtered_measurements = []
        for m in self.measurements:
            if m.alpha == alpha and m.beta == beta and m.sigma_G == sigma_G:
                filtered_measurements += [m]
        return measurement_set(filtered_measurements)
    def get_best_measurement(self,min_chi2=10**6):
        best_measurement = None
        for m in self.measurements:
            if m.total_chi2 < min_chi2:
                min_chi2 = m.total_chi2
                best_measurement = m
        return best_measurement
    def combine_pixels(self):
        z_values = list(set([m.z_value for m in self.measurements]))
        combined_measurements = []
        for z_value in z_values:
            z_set = [m for m in self.measurements if m.z_value == z_value]
            z_parameter_IDs = list(set([m.parameter_ID for m in z_set]))
            for parameter_ID in z_parameter_IDs:
                z_parameter_set = [m for m in z_set if m.parameter_ID == parameter_ID]
                c_m = z_parameter_set[0]
                if len(z_parameter_set) > 1:
                    for m in z_parameter_set[1:]:
                        c_m = measurement.combine_measurements(c_m,m)
                combined_measurements += [c_m]
        return measurement_set(combined_measurements)
    def optimize_s_parameters(self,plot_optimal=False):

        best_measurements = []

        z_values = list(set([m.z_value for m in self.measurements]))
        s_parameter_values_list = list(set([(m.n,m.k1) for m in self.measurements]))

        min_chi2 = 10**6

        for s_parameter_values in s_parameter_values_list:
            n = s_parameter_values[0]
            k1 = s_parameter_values[1]
            s_set = self.s_filter(n,k1)
            print('number measurements in s filtered is:',len(s_set.measurements))
            total_chi2 = 0
            fixed_s_best_measurements = []
            for z_value in z_values:
                z_s_set = s_set.z_filter(z_value)
                best = z_s_set.get_best_measurement()
                total_chi2 += best.total_chi2
                fixed_s_best_measurements += [best]
            if total_chi2 < min_chi2:
                best_measurements = fixed_s_best_measurements
                min_chi2 = total_chi2
        s_optimized_set = measurement_set(measurements=best_measurements)
        print('number measurements in optimised is:',len(s_optimized_set.measurements))

        if plot_optimal:
            plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
            for m in s_optimized_set.measurements:
                plt.errorbar(m.k_kms,m.Pk_kms,fmt='o',label='measured, z={}'.format(m.z_value))
                plt.plot(m.k_kms,P1D_z_kms_PD2013(m.k_kms,m.z_value),label='theory, z={}'.format(m.z_value))
            plt.semilogy()
            plt.semilogx()
            plt.ylabel('Pk1D')
            plt.xlabel('k / kms-1')
            plt.legend()
            plt.grid()
            #plt.savefig('Pk1D_abs_slope1.5_alpha{:2.2f}_beta{:2.2f}_sG{:2.2f}.pdf'.format(alpha,beta,sigma_G_required))
            plt.show()

            plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
            mean_F = []
            alpha = []
            beta = []
            sigma_G = []
            z = []
            for m in s_optimized_set.measurements:
                mean_F += [m.mean_F]
                alpha += [m.alpha]
                beta += [m.beta]
                sigma_G += [m.sigma_G]
                z += [m.z_value]
            plt.errorbar(z,mean_F,fmt='o',label='measured')
            theory_z = np.linspace(1.8,4.0,100)
            plt.plot(theory_z,get_mean_F_model(theory_z),label='theory')
            plt.ylabel('mean F')
            plt.xlabel('z')
            plt.legend()
            plt.grid()
            plt.show()

            plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
            plt.plot(z,alpha,marker='o',label='alpha')
            plt.plot(z,beta,marker='o',label='beta')
            plt.plot(z,sigma_G,marker='o',label='sigma_G')
            plt.ylabel('parameters')
            plt.xlabel('z')
            plt.legend()
            plt.grid()
            plt.show()


        return s_optimized_set



################################################################################
"""
Below: old tuning, no-RSD, theoretical based Method
"""

#Function to return the P1D from Palanque-Delabrouille et al. (2013)
#copied from lyaforecast
def P1D_z_kms_PD2013(k_kms,z,A_F=0.064,B_F=3.55):
    """Fitting formula for 1D P(z,k) from Palanque-Delabrouille et al. (2013).
        Wavenumbers and power in units of km/s. Corrected to be flat at low-k"""
    # numbers from Palanque-Delabrouille (2013)
    #A_F = 0.064
    n_F = -2.55
    alpha_F = -0.1
    #B_F = 3.55
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
