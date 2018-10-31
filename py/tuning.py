import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob

from . import convert, Pk1D, utils

lya = utils.lya_rest

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
    def __init__(self,parameter_ID,z_value,z_width,N_skewers,n,k1,alpha,beta,sigma_G,pixels=[],mean_F=None,k_kms=None,Pk_kms=None,sigma_F=None,cf=None):
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
        self.sigma_F = sigma_F
        self.cf = cf
        return
    def get_details(self):
        details = (self.z_value,self.z_width,self.N_skewers,
                self.n,self.k1,self.alpha,self.beta,self.sigma_G,self.pixels)
        return details
    def add_Pk1D_measurement(self,pixel_object):
        F = pixel_object.lya_absorber.transmission()
        mean_F = np.average(F)
        print('mean F in measuring Pk1D:',mean_F)
        delta_F = F/mean_F - 1
        IVAR = pixel_object.IVAR_rows
        R_hMpc = pixel_object.R
        z = pixel_object.Z
        k_kms, Pk_kms, var_kms = Pk1D.get_Pk1D(delta_F,IVAR,R_hMpc,z,z_value=self.z_value,z_width=self.z_width)
        self.k_kms = k_kms
        self.Pk_kms = Pk_kms
        return
    def add_mean_F_measurement(self,pixel_object):
        self.mean_F = pixel_object.get_mean_flux(pixel_object.lya_absorber,z_value=self.z_value,z_width=self.z_width)
        return
    def add_sigma_F_measurement(self,pixel_object):
        sF = pixel_object.get_sigma_dF(pixel_object.lya_absorber,z_value=self.z_value,z_width=self.z_width)

        Om = 0.3147
        l_hMpc = 0.25
        E_z = np.sqrt(Om*(1+self.z_value)**3 + (1-Om))
        dkms_dhMpc = 100. * E_z / (1+self.z_value)

        # transform to h/Mpc
        k_hMpc = self.k_kms * dkms_dhMpc
        Pk_hMpc = self.Pk_kms / dkms_dhMpc

        # compute Fourier transform of Top-Hat filter of size l_hMpc
        #W_hMpc = np.sinc((k_hMpc*l_hMpc)/(2*np.pi))

        self.sigma_F = np.sqrt((1/np.pi)*np.trapz(Pk_hMpc,k_hMpc)) #(W_hMpc**2)*
        #print('cells: {:2.4f}, hMpc: {:2.4f}, kms: {:2.4f}'.format(sF,self.sigma_F,np.sqrt((1/np.pi)*np.trapz((W_hMpc**2)*self.Pk_kms,self.k_kms))))
        return
    def add_Pk1D_chi2(self,min_k=None,max_k=None,denom="krange10"):
        model_Pk_kms = P1D_z_kms_PD2013(self.z_value,self.k_kms)
        if min_k:
            min_j = max(np.searchsorted(self.k_kms,min_k) - 1,0)
        else:
            min_j = 0
            min_k = 0.
        if max_k:
            max_j = np.searchsorted(self.k_kms,max_k)
        else:
            max_j = -1
            max_k = self.k_kms[-1]
        A = 10**6
        if denom == "uniform5":
            denom = (0.05 * model_Pk_kms)**2
        elif denom == "uniform10":
            denom = (0.10 * model_Pk_kms)**2
        elif denom == "krange5":
            eps = A * np.ones(self.k_kms.shape)
            eps[min_j:max_j] *= 0.05 / A
            denom = (eps * model_Pk_kms)**2
        elif denom == "krange10":
            eps = A * np.ones(self.k_kms.shape)
            eps[min_j:max_j] *= 0.1 / A
            denom = (eps * model_Pk_kms)**2
        elif denom == "krange10_smooth":
            eps = A * np.ones(self.k_kms.shape)
            eps[min_j:max_j] *= 0.1 / A
            smooth_width = (max_k - min_k)/0.01
            lower_smooth = A + (0.1 - A)*np.exp(-((self.k_kms - min_k)**2)/(2*smooth_width**2))
            upper_smooth = A + (0.1 - A)*np.exp(-((self.k_kms - max_k)**2)/(2*smooth_width**2))
            eps[:min_j] = lower_smooth[:min_j]
            eps[max_j:] = upper_smooth[max_j:]
            denom = (eps * model_Pk_kms)**2
        elif denom == "npower":
            k0 = 0.005
            n = 2.
            eps = 0.1 * ((1 + (self.k_kms/k0)**n))
            denom = (eps * model_Pk_kms)**2
        self.Pk_kms_chi2_eps = eps
        chi2 = np.sum(((self.Pk_kms - model_Pk_kms)**2)/denom)
        self.Pk_kms_chi2 = chi2
        return
    def add_mean_F_chi2(self,min_k=None,max_k=None,eps=0.1,mean_F_model='Becker13'):
        model_mean_F = get_mean_F_model(self.z_value,model=mean_F_model)
        denom = (eps * model_mean_F)**2
        chi2 = np.sum(((self.mean_F - model_mean_F)**2)/denom)
        self.mean_F_chi2 = chi2
        return
    def add_sigma_F_chi2(self,min_k=None,max_k=None,eps=0.1,l_hMpc=0.25):
        model_sigma_F = get_sigma_dF_P1D(self.z_value,l_hMpc=l_hMpc)
        denom = (eps * model_sigma_F)**2
        chi2 = np.sum(((self.sigma_F - model_sigma_F)**2)/denom)
        self.sigma_F_chi2 = chi2
        #print(self.z_value,model_sigma_F,self.sigma_F)
        return
    def add_total_chi2(self):
        chi2 = self.Pk_kms_chi2 + self.mean_F_chi2
        self.total_chi2 = chi2
        return
    @classmethod
    def combine_measurements(cls,m1,m2):
        if utils.confirm_identical(m1.parameter_ID,m2.parameter_ID,item_name='parameter_ID'):
            parameter_ID = m1.parameter_ID
            n = m1.n
            k1 = m1.k1
            alpha = m1.alpha
            beta = m1.beta
            sigma_G = m1.sigma_G
        if utils.confirm_identical(m1.z_value,m2.z_value,item_name='z_value'):
            z_value = m1.z_value
        if utils.confirm_identical(m1.z_width,m2.z_width,item_name='z_width'):
            z_width = m1.z_width
        N_skewers = m1.N_skewers + m2.N_skewers
        pixels = m1.pixels + m2.pixels
        mean_F = (m1.mean_F*m1.N_skewers + m2.mean_F*m2.N_skewers)/(m1.N_skewers + m2.N_skewers)
        k_kms_check = True
        if utils.confirm_identical(m1.k_kms,m2.k_kms,item_name='k_kms',array=True):
            k_kms = m1.k_kms
        Pk_kms = (m1.Pk_kms*m1.N_skewers + m2.Pk_kms*m2.N_skewers)/(m1.N_skewers + m2.N_skewers)
        sigma_F = np.sqrt(((m1.sigma_F**2)*m1.N_skewers + (m2.sigma_F**2)*m2.N_skewers)/(m1.N_skewers + m2.N_skewers))
        #May need to work on this?
        cf = None
        return measurement(parameter_ID,z_value,z_width,N_skewers,n,k1,alpha,beta,sigma_G,pixels=pixels,mean_F=mean_F,k_kms=k_kms,Pk_kms=Pk_kms,sigma_F=sigma_F,cf=cf)
    @classmethod
    def load_measurement(cls,hdu):
        parameter_ID = hdu.header['param_ID']
        z_value = hdu.header['z_value']
        z_width = hdu.header['z_width']
        N_skewers = hdu.header['N_skw']

        n = hdu.header['n']
        k1 = hdu.header['k1']
        alpha = hdu.header['alpha']
        beta = hdu.header['beta']
        sigma_G = hdu.header['sigma_G']

        #Not sure how to do this...
        N_pix = hdu.header['N_pix']
        if hdu.header['pixels'] == 'multiple':
            pixels = ['']*N_pix
        else:
            pixels = [hdu.header['pixels']]

        mean_F = hdu.header['mean_F']
        k_kms = hdu.data['k_kms']
        Pk_kms = hdu.data['Pk_kms']
        cf = hdu.header['cf']
        return cls(parameter_ID,z_value,z_width,N_skewers,n,k1,alpha,beta,sigma_G,pixels=pixels,mean_F=mean_F,k_kms=k_kms,Pk_kms=Pk_kms,cf=cf)
    def make_HDU(self):
        header = fits.Header()
        header['param_ID'] = self.parameter_ID
        header['z_value'] = self.z_value
        header['z_width'] = self.z_width
        header['N_skw'] = self.N_skewers
        header['n'] = self.n
        header['k1'] = self.k1
        header['alpha'] = self.alpha
        header['beta'] = self.beta
        header['sigma_G'] = self.sigma_G
        header['N_pix'] = len(self.pixels)
        if len(self.pixels) > 1:
            header['pixels'] = 'multiple'
        else:
            header['pixels'] = self.pixels[0]
        header['mean_F'] = self.mean_F
        header['cf'] = self.cf

        Pk_data = list(zip(self.k_kms,self.Pk_kms))
        dtype = [('k_kms', 'f8'), ('Pk_kms', 'f8')]
        measurement_1 = np.array(Pk_data,dtype=dtype)
        cols = fits.ColDefs(measurement_1)
        hdu = fits.BinTableHDU.from_columns(cols,header=header,name='Pk1D')
        return hdu
    # TODO: functions to make plots

class function_measurement:
    def __init__(self,parameter_ID,z_value,z_width,N_skewers,n,k1,C0,C1,C2,beta,D0,D1,D2,pixels=[],mean_F=None,k_kms=None,Pk_kms=None,sigma_F=None,cf=None):
        self.parameter_ID = parameter_ID
        self.z_value = z_value
        self.z_width = z_width
        self.N_skewers = N_skewers

        self.n = n
        self.k1 = k1
        self.C0 = C0
        self.C1 = C1
        self.C2 = C2
        self.beta = beta
        self.D0 = D0
        self.D1 = D1
        self.D2 = D2

        self.pixels = pixels

        self.mean_F = mean_F
        self.k_kms = k_kms
        self.Pk_kms = Pk_kms
        self.sigma_F = sigma_F
        self.cf = cf
        return
    def get_details(self):
        details = (self.z_value,self.z_width,self.N_skewers,
                self.n,self.k1,self.C0,self.C1,self.C2,self.D0,self.D1,self.D2,self.beta,self.pixels)
        return details
    def add_Pk1D_measurement(self,pixel_object):
        if not self.mean_F:
            self.add_mean_F_measurement(pixel_object)
        mean_F = self.mean_F
        F = pixel_object.lya_absorber.transmission()
        delta_F = F/mean_F - 1
        IVAR = pixel_object.IVAR_rows
        R_hMpc = pixel_object.R
        z = pixel_object.Z
        k_kms, Pk_kms, var_kms = Pk1D.get_Pk1D(delta_F,IVAR,R_hMpc,z,z_value=self.z_value,z_width=self.z_width)
        self.k_kms = k_kms
        self.Pk_kms = Pk_kms
        return
    def add_mean_F_measurement(self,pixel_object):
        
        self.mean_F = pixel_object.get_mean_flux(pixel_object.lya_absorber,z_value=self.z_value,z_width=self.z_width)
        
        return
    def add_sigma_F_measurement(self,pixel_object):
        sF = pixel_object.get_sigma_dF(pixel_object.lya_absorber,z_value=self.z_value,z_width=self.z_width)

        Om = 0.3147
        l_hMpc = 0.25
        E_z = np.sqrt(Om*(1+self.z_value)**3 + (1-Om))
        dkms_dhMpc = 100. * E_z / (1+self.z_value)

        # transform to h/Mpc
        k_hMpc = self.k_kms * dkms_dhMpc
        Pk_hMpc = self.Pk_kms / dkms_dhMpc

        # compute Fourier transform of Top-Hat filter of size l_hMpc
        #W_hMpc = np.sinc((k_hMpc*l_hMpc)/(2*np.pi))

        self.sigma_F = np.sqrt((1/np.pi)*np.trapz(Pk_hMpc,k_hMpc)) #(W_hMpc**2)*
        #print('cells: {:2.4f}, hMpc: {:2.4f}, kms: {:2.4f}'.format(sF,self.sigma_F,np.sqrt((1/np.pi)*np.trapz((W_hMpc**2)*self.Pk_kms,self.k_kms))))
        return
    def add_Pk1D_chi2(self,min_k=None,max_k=None,denom="krange10"):
        model_Pk_kms = P1D_z_kms_PD2013(self.z_value,self.k_kms)
        if min_k:
            min_j = max(np.searchsorted(self.k_kms,min_k) - 1,0)
        else:
            min_j = 0
            min_k = 0.
        if max_k:
            max_j = np.searchsorted(self.k_kms,max_k)
        else:
            max_j = -1
            max_k = self.k_kms[-1]
        A = 10**6
        if denom == "uniform5":
            denom = (0.05 * model_Pk_kms)**2
        elif denom == "uniform10":
            denom = (0.10 * model_Pk_kms)**2
        elif denom == "krange5":
            eps = A * np.ones(self.k_kms.shape)
            eps[min_j:max_j] *= 0.05 / A
            denom = (eps * model_Pk_kms)**2
        elif denom == "krange10":
            eps = A * np.ones(self.k_kms.shape)
            eps[min_j:max_j] *= 0.1 / A
            denom = (eps * model_Pk_kms)**2
        elif denom == "krange10_smooth":
            eps = A * np.ones(self.k_kms.shape)
            eps[min_j:max_j] *= 0.1 / A
            smooth_width = (max_k - min_k)/0.01
            lower_smooth = A + (0.1 - A)*np.exp(-((self.k_kms - min_k)**2)/(2*smooth_width**2))
            upper_smooth = A + (0.1 - A)*np.exp(-((self.k_kms - max_k)**2)/(2*smooth_width**2))
            eps[:min_j] = lower_smooth[:min_j]
            eps[max_j:] = upper_smooth[max_j:]
            denom = (eps * model_Pk_kms)**2
        elif denom == "npower":
            k0 = 0.005
            n = 2.
            eps = 0.1 * ((1 + (self.k_kms/k0)**n))
            denom = (eps * model_Pk_kms)**2
        self.Pk_kms_chi2_eps = eps
        chi2 = np.sum(((self.Pk_kms - model_Pk_kms)**2)/denom)
        self.Pk_kms_chi2 = chi2
        return
    def add_mean_F_chi2(self,min_k=None,max_k=None,eps=0.1,mean_F_model='Becker13'):
        model_mean_F = get_mean_F_model(self.z_value,model=mean_F_model)
        denom = (eps * model_mean_F)**2
        chi2 = np.sum(((self.mean_F - model_mean_F)**2)/denom)
        self.mean_F_chi2 = chi2
        return
    def add_sigma_F_chi2(self,min_k=None,max_k=None,eps=0.1,l_hMpc=0.25):
        model_sigma_F = get_sigma_dF_P1D(self.z_value,l_hMpc=l_hMpc)
        denom = (eps * model_sigma_F)**2
        chi2 = np.sum(((self.sigma_F - model_sigma_F)**2)/denom)
        self.sigma_F_chi2 = chi2
        #print(self.z_value,model_sigma_F,self.sigma_F)
        return
    def add_total_chi2(self):
        chi2 = self.Pk_kms_chi2 + self.mean_F_chi2
        self.total_chi2 = chi2
        return
    @classmethod
    def combine_measurements(cls,m1,m2):
        if utils.confirm_identical(m1.parameter_ID,m2.parameter_ID,item_name='parameter_ID'):
            parameter_ID = m1.parameter_ID
            n = m1.n
            k1 = m1.k1
            C0 = m1.C0
            C1 = m1.C1
            C2 = m1.C2
            beta = m1.beta
            D0 = m1.D0
            D1 = m1.D1
            D2 = m1.D2
        if utils.confirm_identical(m1.z_value,m2.z_value,item_name='z_value'):
            z_value = m1.z_value
        if utils.confirm_identical(m1.z_width,m2.z_width,item_name='z_width'):
            z_width = m1.z_width
        N_skewers = m1.N_skewers + m2.N_skewers
        pixels = m1.pixels + m2.pixels
        mean_F = (m1.mean_F*m1.N_skewers + m2.mean_F*m2.N_skewers)/(m1.N_skewers + m2.N_skewers)
        k_kms_check = True
        if utils.confirm_identical(m1.k_kms,m2.k_kms,item_name='k_kms',array=True):
            k_kms = m1.k_kms
        Pk_kms = (m1.Pk_kms*m1.N_skewers + m2.Pk_kms*m2.N_skewers)/(m1.N_skewers + m2.N_skewers)
        sigma_F = np.sqrt(((m1.sigma_F**2)*m1.N_skewers + (m2.sigma_F**2)*m2.N_skewers)/(m1.N_skewers + m2.N_skewers))
        #May need to work on this?
        cf = None
        return function_measurement(parameter_ID,z_value,z_width,N_skewers,n,k1,C0,C1,C2,beta,D0,D1,D2,pixels=pixels,mean_F=mean_F,k_kms=k_kms,Pk_kms=Pk_kms,sigma_F=sigma_F,cf=cf)
    @classmethod
    def load_measurement(cls,hdu):
        parameter_ID = hdu.header['param_ID']
        z_value = hdu.header['z_value']
        z_width = hdu.header['z_width']
        N_skewers = hdu.header['N_skw']

        n = hdu.header['n']
        k1 = hdu.header['k1']
        alpha = hdu.header['alpha']
        beta = hdu.header['beta']
        sigma_G = hdu.header['sigma_G']

        #Not sure how to do this...
        N_pix = hdu.header['N_pix']
        if hdu.header['pixels'] == 'multiple':
            pixels = ['']*N_pix
        else:
            pixels = [hdu.header['pixels']]

        mean_F = hdu.header['mean_F']
        k_kms = hdu.data['k_kms']
        Pk_kms = hdu.data['Pk_kms']
        cf = hdu.header['cf']
        return cls(parameter_ID,z_value,z_width,N_skewers,n,k1,alpha,beta,sigma_G,pixels=pixels,mean_F=mean_F,k_kms=k_kms,Pk_kms=Pk_kms,cf=cf)
    def make_HDU(self):
        header = fits.Header()
        header['param_ID'] = self.parameter_ID
        header['z_value'] = self.z_value
        header['z_width'] = self.z_width
        header['N_skw'] = self.N_skewers
        header['n'] = self.n
        header['k1'] = self.k1
        header['alpha'] = self.alpha
        header['beta'] = self.beta
        header['sigma_G'] = self.sigma_G
        header['N_pix'] = len(self.pixels)
        if len(self.pixels) > 1:
            header['pixels'] = 'multiple'
        else:
            header['pixels'] = self.pixels[0]
        header['mean_F'] = self.mean_F
        header['cf'] = self.cf

        Pk_data = list(zip(self.k_kms,self.Pk_kms))
        dtype = [('k_kms', 'f8'), ('Pk_kms', 'f8')]
        measurement_1 = np.array(Pk_data,dtype=dtype)
        cols = fits.ColDefs(measurement_1)
        hdu = fits.BinTableHDU.from_columns(cols,header=header,name='Pk1D')
        return hdu
    # TODO: functions to make plots

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
                        c_m = function_measurement.combine_measurements(c_m,m)
                combined_measurements += [c_m]
        return measurement_set(combined_measurements)
    def optimize_s_parameters(self,plot_optimal=False,mean_F_model='Becker13'):

        best_measurements = []

        z_values = list(set([m.z_value for m in self.measurements]))
        s_parameter_values_list = list(set([(m.n,m.k1) for m in self.measurements]))

        min_chi2 = 10**6

        for s_parameter_values in s_parameter_values_list:
            n = s_parameter_values[0]
            k1 = s_parameter_values[1]
            print(n,k1)
            s_set = self.s_filter(n,k1)
            print('number measurements in s filtered is:',len(s_set.measurements))
            total_chi2 = 0
            fixed_s_best_measurements = []
            for z_value in z_values:
                print('->',z_value)
                z_s_set = s_set.z_filter(z_value)
                best = z_s_set.get_best_measurement()
                print('->-> chi2 Pk1D {:2.2f}, mean_F {:2.2f}, total {:2.2f}'.format(best.Pk_kms_chi2,best.mean_F_chi2,best.total_chi2))
                print('->-> alpha {:2.2f}, beta {:2.2f}, sigma_G {:2.2f}'.format(best.alpha,best.beta,best.sigma_G))
                total_chi2 += best.total_chi2
                fixed_s_best_measurements += [best]
            if total_chi2 < min_chi2:
                best_measurements = fixed_s_best_measurements
                min_chi2 = total_chi2
            print(' ')
        s_optimized_set = measurement_set(measurements=best_measurements)
        print('number measurements in optimised is:',len(s_optimized_set.measurements))
        for best in s_optimized_set.measurements:
            print(best.z_value)
            print('->-> chi2 Pk1D {:2.2f}, mean_F {:2.2f}, total {:2.2f}'.format(best.Pk_kms_chi2,best.mean_F_chi2,best.total_chi2))
            print('->-> alpha {:2.2f}, beta {:2.2f}, sigma_G {:2.2f}'.format(best.alpha,best.beta,best.sigma_G))
            print('->-> n {:2.2f}, k1 {:2.6f}'.format(best.n,best.k1))
        if plot_optimal:
            plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
            for m in s_optimized_set.measurements:
                plt.errorbar(m.k_kms,m.Pk_kms,fmt='o',label='measured, z={}'.format(m.z_value))
                plt.plot(m.k_kms,P1D_z_kms_PD2013(m.z_value,m.k_kms),label='theory, z={}'.format(m.z_value))
            plt.semilogy()
            plt.semilogx()
            plt.ylabel('Pk1D')
            plt.xlabel('k / kms-1')
            plt.legend()
            plt.grid()
            plt.savefig('Pk1D_010.pdf')
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
            plt.plot(theory_z,get_mean_F_model(theory_z,model=mean_F_model),label='theory')
            plt.ylabel('mean F')
            plt.xlabel('z')
            plt.legend()
            plt.grid()
            plt.savefig('mean_F_010.pdf')
            plt.show()

            plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
            plt.plot(z,alpha,marker='o',label='alpha')
            plt.plot(z,beta,marker='o',label='beta')
            plt.plot(z,sigma_G,marker='o',label='sigma_G')
            plt.ylabel('parameters')
            plt.xlabel('z')
            plt.legend()
            plt.grid()
            plt.savefig('parameters_010.pdf')
            plt.show()

        n_grids = len(z_values) * 3
        n_values = np.sort(list(set([spv[0] for spv in s_parameter_values_list])))
        k1_values = np.sort(list(set([spv[1] for spv in s_parameter_values_list])))
        colour_grids = np.zeros((n_grids,len(n_values),len(k1_values)))
        for s_parameter_values in s_parameter_values_list:
            n = s_parameter_values[0]
            k1 = s_parameter_values[1]
            n_i = np.searchsorted(n_values,n)
            k1_i = np.searchsorted(k1_values,k1)
            s_set = self.s_filter(n,k1)
            total_chi2 = 0
            for j,z_value in enumerate(z_values):
                z_s_set = s_set.z_filter(z_value)
                best = z_s_set.get_best_measurement()
                colour_grids[j*3+0,k1_i,n_i] = best.Pk_kms_chi2
                colour_grids[j*3+1,k1_i,n_i] = best.mean_F_chi2
                colour_grids[j*3+2,k1_i,n_i] = best.total_chi2

        for k,z_value in enumerate(z_values):
            fig, ax = plt.subplots(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
            im = ax.imshow(colour_grids[k*3+0,:,:],cmap='YlGn',vmin=0,vmax=np.max(colour_grids))
            # Create colorbar
            cbar = ax.figure.colorbar(im, ax=ax)
            cbar.ax.set_ylabel("chi2", rotation=-90, va="bottom")
            # We want to show all ticks...
            ax.set_xticks(np.arange(len(n_values)))
            ax.set_yticks(np.arange(len(k1_values)))
            # ... and label them with the respective list entries
            ax.set_xticklabels(n_values)
            ax.set_yticklabels(k1_values)
            # Rotate the tick labels and set their alignment.
            plt.setp(ax.get_xticklabels(), ha="right")
            # Loop over data dimensions and create text annotations.
            for i in range(len(n_values)):
                for j in range(len(k1_values)):
                    text = ax.text(j, i, round(colour_grids[k*3+0,i,j],2),
                                   ha="center", va="center", color=(0,0,0))
            ax.set_title('Pk chi2 values, z={}'.format(z_value))
            plt.savefig('colour_z{}_{}.pdf'.format(z_value,'Pk'))
            plt.show()

            fig, ax = plt.subplots(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
            im = ax.imshow(colour_grids[k*3+1,:,:],cmap='YlGn',vmin=0,vmax=np.max(colour_grids))
            # Create colorbar
            cbar = ax.figure.colorbar(im, ax=ax)
            cbar.ax.set_ylabel("chi2", rotation=-90, va="bottom")
            # We want to show all ticks...
            ax.set_xticks(np.arange(len(n_values)))
            ax.set_yticks(np.arange(len(k1_values)))
            # ... and label them with the respective list entries
            ax.set_xticklabels(n_values)
            ax.set_yticklabels(k1_values)
            # Rotate the tick labels and set their alignment.
            plt.setp(ax.get_xticklabels(), ha="right")
            # Loop over data dimensions and create text annotations.
            for i in range(len(n_values)):
                for j in range(len(k1_values)):
                    text = ax.text(j, i, round(colour_grids[k*3+1,i,j],2),
                                   ha="center", va="center", color=(0,0,0))
            ax.set_title('mean F chi2 values, z={}'.format(z_value))
            plt.savefig('colour_z{}_{}.pdf'.format(z_value,'mean_F'))
            plt.show()

            fig, ax = plt.subplots(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
            im = ax.imshow(colour_grids[k*3+2,:,:],cmap='YlGn',vmin=0,vmax=np.max(colour_grids))
            # Create colorbar
            cbar = ax.figure.colorbar(im, ax=ax)
            cbar.ax.set_ylabel("chi2", rotation=-90, va="bottom")
            # We want to show all ticks...
            ax.set_xticks(np.arange(len(n_values)))
            ax.set_yticks(np.arange(len(k1_values)))
            # ... and label them with the respective list entries
            ax.set_xticklabels(n_values)
            ax.set_yticklabels(k1_values)
            # Rotate the tick labels and set their alignment.
            plt.setp(ax.get_xticklabels(), ha="right")
            # Loop over data dimensions and create text annotations.
            for i in range(len(n_values)):
                for j in range(len(k1_values)):
                    text = ax.text(j, i, round(colour_grids[k*3+2,i,j],2),
                                   ha="center", va="center", color=(0,0,0))
            ax.set_title('Total chi2 values, z={}'.format(z_value))
            plt.savefig('colour_z{}_{}.pdf'.format(z_value,'total'))
            plt.show()

        if len(z_values) > 1:
            collapsed_z_colour_grids = np.zeros((3,len(n_values),len(k1_values)))
            for j,z_value in enumerate(z_values):
                collapsed_z_colour_grids[0,:,:] += colour_grids[3*j+0,:,:]
                collapsed_z_colour_grids[1,:,:] += colour_grids[3*j+1,:,:]
                collapsed_z_colour_grids[2,:,:] += colour_grids[3*j+2,:,:]

            fig, ax = plt.subplots(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
            im = ax.imshow(collapsed_z_colour_grids[0,:,:],cmap='YlGn',vmin=0,vmax=np.max(colour_grids))
            # Create colorbar
            cbar = ax.figure.colorbar(im, ax=ax)
            cbar.ax.set_ylabel("chi2", rotation=-90, va="bottom")
            # We want to show all ticks...
            ax.set_xticks(np.arange(len(n_values)))
            ax.set_yticks(np.arange(len(k1_values)))
            # ... and label them with the respective list entries
            ax.set_xticklabels(n_values)
            ax.set_yticklabels(k1_values)
            # Rotate the tick labels and set their alignment.
            plt.setp(ax.get_xticklabels(), ha="right")
            # Loop over data dimensions and create text annotations.
            for i in range(len(n_values)):
                for j in range(len(k1_values)):
                    text = ax.text(j, i, round(collapsed_z_colour_grids[0,i,j],2),
                                   ha="center", va="center", color=(0,0,0))
            ax.set_title('Pk chi2 values, all z values')
            plt.savefig('colour_{}.pdf'.format('Pk'))
            plt.show()

            fig, ax = plt.subplots(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
            im = ax.imshow(collapsed_z_colour_grids[1,:,:],cmap='YlGn',vmin=0,vmax=np.max(colour_grids))
            # Create colorbar
            cbar = ax.figure.colorbar(im, ax=ax)
            cbar.ax.set_ylabel("chi2", rotation=-90, va="bottom")
            # We want to show all ticks...
            ax.set_xticks(np.arange(len(n_values)))
            ax.set_yticks(np.arange(len(k1_values)))
            # ... and label them with the respective list entries
            ax.set_xticklabels(n_values)
            ax.set_yticklabels(k1_values)
            # Rotate the tick labels and set their alignment.
            plt.setp(ax.get_xticklabels(), ha="right")
            # Loop over data dimensions and create text annotations.
            for i in range(len(n_values)):
                for j in range(len(k1_values)):
                    text = ax.text(j, i, round(collapsed_z_colour_grids[1,i,j],2),
                                   ha="center", va="center", color=(0,0,0))
            ax.set_title('mean F chi2 values, all z values')
            plt.savefig('colour_{}.pdf'.format('mean_F'))
            plt.show()

            fig, ax = plt.subplots(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
            im = ax.imshow(collapsed_z_colour_grids[2,:,:],cmap='YlGn',vmin=0,vmax=np.max(colour_grids))
            # Create colorbar
            cbar = ax.figure.colorbar(im, ax=ax)
            cbar.ax.set_ylabel("chi2", rotation=-90, va="bottom")
            # We want to show all ticks...
            ax.set_xticks(np.arange(len(n_values)))
            ax.set_yticks(np.arange(len(k1_values)))
            # ... and label them with the respective list entries
            ax.set_xticklabels(n_values)
            ax.set_yticklabels(k1_values)
            # Rotate the tick labels and set their alignment.
            plt.setp(ax.get_xticklabels(), ha="right")
            # Loop over data dimensions and create text annotations.
            for i in range(len(n_values)):
                for j in range(len(k1_values)):
                    text = ax.text(j, i, round(collapsed_z_colour_grids[2,i,j],2),
                                   ha="center", va="center", color=(0,0,0))
            ax.set_title('Total chi2 values, all z values')
            plt.savefig('colour_{}.pdf'.format('total'))
            plt.show()

        return s_optimized_set
    def save(self,filepath,existing='overwrite'):
        #Make a list of HDUs from the current set of measurements.
        list_of_hdus = []
        for m in self.measurements:
            list_of_hdus += [m.make_HDU()]

        #Check if a file already exists. If so, follow the chosen path, if not, make a new file.
        if len(glob.glob(filepath)) > 0:
            if existing == 'combine':
                #Open the existing file and get the list of what parameter sets it includes.
                hdulist = fits.open(filepath)
                existing_parameters_list = []
                for existing_hdu in hdulist[1:]:
                    m = measurement.load_measurement(existing_hdu)
                    existing_parameters_list += [m.get_details()]
                    print('#####')
                    print(existing_parameters_list)

                for i,hdu in enumerate(list_of_hdus):
                    details = measurement.load_measurement(hdu).get_details()
                    print(details)
                    if details in existing_parameters_list:
                        print('measurement with parameters {} already included in {}, existing measurement retained'.format(details,filepath))
                    else:
                        hdulist.append(hdu)
            elif existing == 'overwrite':
                prihdr = fits.Header()
                prihdu = fits.PrimaryHDU(header=prihdr)
                hdulist = fits.HDUList([prihdu]+list_of_hdus)
            else:
                raise ValueError('File already exists, specified behaviour not recognised')
        else:
            prihdr = fits.Header()
            prihdu = fits.PrimaryHDU(header=prihdr)
            hdulist = fits.HDUList([prihdu]+list_of_hdus)
        hdulist.writeto(filepath,overwrite=True)
        hdulist.close()
        return

def fit_function_to_data(x,y,new_x):

    from iminuit import Minuit

    def fitting_function(x,A0,A1,A2):
        return A0 * (x ** A1) + A2

    def least_squares(A0,A1,A2):
        return np.sum((y - fitting_function(x,A0,A1,A2))**2)

    m = Minuit(least_squares,A0=1.,A1=1.,A2=1.,error_A0=0.1,error_A1=0.1,error_A2=0.1,errordef=1)
    fmin,param = m.migrad()

    new_y = fitting_function(new_x,**m.values)

    return new_y

################################################################################
"""
Below: old tuning, no-RSD, theoretical based Method
"""

#Function to return the P1D from Palanque-Delabrouille et al. (2013)
#copied from lyaforecast
def P1D_z_kms_PD2013(z,k_kms,A_F=0.064,B_F=3.55):
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
def get_mean_F_model(z,model='Becker13'):
    if model == 'Becker13':
        mean_F = np.exp(-0.751*(((1+z)/4.5)**2.9)+0.132)
    elif model == 'FontRibera12':
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
def get_sigma_dF_P1D(z,l_hMpc=0.25,Om=0.3147):
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
