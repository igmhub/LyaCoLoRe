import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import json
from scipy.stats import truncnorm

from lyacolore import bias, convert, Pk1D, transformation, utils

def get_initial_transformation_and_minuit_input(filepath,random=False,seed=0):

    t = transformation.Transformation.make_transformation_from_tuning_file(filepath,random=False,seed=0)
    minuit_input_dict, param_names = tuning_file_to_minuit_input(filepath,random=False,seed=0)

    return t, minuit_input_dict, param_names

def tuning_file_to_minuit_input(filepath,random=False,seed=0):

    # Open the tuning file.
    with open(filepath, 'r') as json_file:
        data = json.load(json_file)

    # Randomise values if desired.
    if random:
        transformation.randomise_parameter_values(data,seed=seed)

    minuit_input_dict = {}
    param_names = []

    # For each tuning quantity:
    for k in data.keys():
        # For each parameter defining that quantity:
        for p in data[k]['parameters'].keys():
            p_dict = data[k]['parameters'][p]

            # Make a parameter name.
            p_name = k+'-'+p
            param_names += [p_name]

            # Add the parameter's initial value, error and fix to the minuit
            # input dictionary.
            minuit_input_dict[p_name] = p_dict['value']
            minuit_input_dict['error_'+p_name] = p_dict['error']
            minuit_input_dict['fix_'+p_name] = p_dict['fix']

            # If desired, also add limits.
            if ('lower_limit' in p_dict.keys()) and ('upper_limit' in p_dict.keys()):
                limit = (p_dict['lower_limit'], p_dict['upper_limit'])
            elif ('lower_limit' in p_dict.keys()):
                limit = (p_dict['lower_limit'], float("infinity"))
            elif ('upper_limit' in p_dict.keys()):
                limit = (-float("infinity"), p_dict['upper_limit'])
            else:
                limit = (-float("infinity"), float("infinity"))
            minuit_input_dict['limit_'+p_name] = limit

    return minuit_input_dict, param_names

class MinimisationObject:
    def __init__(self,tuning_args,run_args):
        ## Given the inputs from config files, this constructs a set of SimulationData objects.

        ## Choose the QSOs from the master file.
        master = fits.open(run_args.out_dir+'/master.fits')
        nobj_total = len(master[1].data)
        w = np.random.choice(range(nobj_total),tuning_args.n_skewers,replace=False)
        input_objects = master[1].data[w]
        master.close()

        ## Load the CoLoRe input skewers.
        fnums = np.sort(list(set(input_objects['FILENUM'])))
        fpaths = {fnum:get_in_file_name(run_args.in_dir,run_args.in_file_prefix,fnum) for fnum in fnums}
        fmockids = {fnum:input_objects['MOCKID'][(input_objects['FILENUM']==fnum)] for filenum in filenums}
        tasks = [fpaths[fnum],fnum,run_args.file_format,run_args.skewer_type,filemockids[fnum],run_args.lambda_min,run_args.rest_frame_weights_cut,None]
        if __name__ == '__main__':
            pool = Pool(processes = tuning_args.nproc)
            results = []
            for task in tasks:
                pool.apply_async(simulation_data.get_skewers_object,task,callback=log_result,error_callback=log_error)
            pool.close()
            pool.join()

        ## Assign the results as an attribute.
        self.data = results

        ## Load the initial tuning file as minuit input.
        minuit_input = tuning_file_to_minuit_input(
            file=tuning_args.initial_parameter_file,
            random=tuning_args.randomise_initial_parameter_values,
            seed=tuning_args.seed,
            )
        self.minuit_input = minuit_input

        ## Check to see if the velocity boost is constant.
        fixed_velocity_boost = True
        for k in minuit_input.keys():
            if (k[-3:]=='a_v') & (k[-4]=='_') & (k[:3]=='fix'):
                fixed_velocity_boost *= minuit_input[k]

        ## If so, pre-calculate the RSDs to make things faster.
        if fixed_velocity_boost:
            self.add_transformation(tuning_args.initial_parameter_file)
            self.calculate_rsds(tuning_args,run_args)

        return

    def add_transformation(self,parameter_file):

        transformation = Transformation.make_transformation_from_file(parameter_file)

        for d in self.data:
            d.add_transformation(transformation)
            d.scale_velocities(use_transformation=True)

        return

    def calculate_rsds(self,tuning_args,run_args):

        tasks = range(len(self.data))

        def add_rsds_to_data_object(data):

            seed = int(pixel * 10**5 + args.seed)

            #trim skewers to the minimal length
            lambda_buffer = 100. #Angstroms
            z_lower_cut = np.min(tuning_args.z_values) - tuning_args.z_width/2.
            z_upper_cut = np.max(tuning_args.z_values) + tuning_args.z_width/2.
            lambda_min_val = np.min([run_args.lambda_min,utils.lya_rest*(1 + z_lower_cut)]) - lambda_buffer
            lambda_max_val = utils.lya_rest*(1 + z_upper_cut) + lambda_buffer
            data.trim_skewers(lambda_min_val,run_args.min_cat_z,lambda_max=lambda_max_val,whole_lambda_range=False)

            #Add small scale fluctuations to the skewers.
            generator = np.random.RandomState(seed)
            data.add_small_scale_fluctuations(args.cell_size,generator,white_noise=False,lambda_min=0.0,IVAR_cutoff=args.lambda_rest_max,use_transformation=True,remove_P1D_data=remove_P1D_data)

            #print('{:3.2f} checkpoint extra power'.format(time.time()-t))
            t = time.time()

            #If needed, compute the physical skewers
            if args.skewer_type == 'gaussian':
                data.compute_physical_skewers()

            #Compute the tau skewers and add RSDs
            data.compute_tau_skewers(data.lya_absorber)
            #print('{:3.2f} checkpoint tau'.format(time.time()-t))
            t = time.time()

            if prep:
                data.compute_RSD_weights(thermal=False)

        return

    def __call__(self,**args):
        ## This calculates a loss function to be minimised, given a set of
        ## parameters.
        return

    def save_tuning_file(self,minuit_output,tuning_filepath):

        tuning_dict = {}
        tuning_quantities = ['tau0', 'alpha', 'seps', 'ssf']

        # For each tuning quantity, put the information for parameters relating to
        # that quantity into a list. Store that list in a {tuning quantity: list}
        # dictionary.
        tuning_pre_dict = {}
        for tq in tuning_quantities:
            tq_list = []
            for parameter_info in minuit_output.get_param_states():
                if (parameter_info['name'][:len(tq)] == tq) and (parameter_info['name'][len(tq):len(tq)+1] == '_'):
                    tq_list += [parameter_info]
            tuning_pre_dict[tq] = tq_list

        # Check that all parameters from the output are accounted for.
        assert np.sum([len(tuning_pre_dict[tq]) for tq in tuning_quantities]) == len(minuit_output.get_param_states())

        # Arrange the tuning output dictionary.


        # Save as a json file.

        return


################################################################################
"""
Below: new tuning, measurement based
"""



"""class TuningQuantity:
    def __init__(self):
        return

    @classmethod
    def make_tuning_quantity_from_dict(cls,tq_dict):
        # Instantiate class.
        tuning_quantity = cls()

        # Load the parameters' properties.
        parameters = {k:tq_dict['parameters'][k]['value'] for k in tq_dict['parameters'].keys()}
        tuning_quantity = functiontypes[tq_dict['functiontype']](**parameters)

        return"""

"""
    #Function to add the transformation functions by interpolating data.
    def add_parameters_from_data(self,z_values,tau0_values,texp_values,seps_values,n,k1,R_kms,a_v):
        # For each of the redshift-dependent parameters, define a function which
        # interpolates the data given and add it to the object.
        def f_tau0_z(z):
            return np.exp(np.interp(np.log(z),np.log(z_values),np.log(tau0_values)))
        def f_texp_z(z):
            return np.exp(np.interp(np.log(z),np.log(z_values),np.log(texp_values)))
        def f_seps_z(z):
            return np.exp(np.interp(np.log(z),np.log(z_values),np.log(seps_values)))
        self.add_zdep_parameters_from_functions(f_tau0_z,f_texp_z,f_seps_z)

        # Add each of the single-value parameters to the object.
        self.add_singval_parameters(n=n,k1=k1,R_kms=R_kms,a_v=a_v)

        return

    # TODO: UPDATE TUNING FILE COLUMN NAMES
    #Function to add the transformation functions by interpolating data from a
    #given tuning file.
    def add_parameters_from_file(self,filepath):

        # Open the file.
        h = fits.open(filepath)

        # Get the redshift-dependent parameter arrays.
        z = h[1].data['z']
        try:
            tau0 = h[1].data['tau0']
        except:
            tau0 = h[1].data['alpha']
        try:
            texp = h[1].data['texp']
        except:
            texp = h[1].data['beta']
        try:
            seps = h[1].data['seps']
        except:
            seps = h[1].data['sigma_G']

        # Get the single-value parameters from the header.
        n = h[1].header['n']
        k1 = h[1].header['k1']
        R_kms = h[1].header['R']
        try:
            a_v = h[1].header['a_v']
        except:
            a_v = h[1].header['vb']

        self.add_parameters_from_data(z,tau0,texp,seps,n,k1,R_kms,a_v)

        return

    #Function to add the z dependent parameter functions by giving the functions.
    def add_zdep_parameters_from_functions(self,f_tau0_z,f_texp_z,f_seps_z):
        self.f_tau0_z = f_tau0_z
        self.f_texp_z = f_texp_z
        self.f_seps_z = f_seps_z
        return

    #Function to add the single value parameters to the object.
    def add_singval_parameters(self,n=0.7,k1=0.001,R_kms=25.,a_v=1.):
        self.n = n
        self.k1 = k1
        self.R_kms = R_kms
        self.a_v = a_v
        return

    #Function to evaluate tau0, the normalisation of the FGPA.
    def get_tau0(self,z):
        return self.f_tau0_z(z)

    #Function to evaluate texp, the exponent of the FGPA.
    def get_texp(self,z):
        return self.f_texp_z(z)

    #Function to evaluate sigma_epsilon, the std of the extra power.
    def get_seps(self,z):
        return self.f_seps_z(z)"""

class function_measurement:
    def __init__(self,parameter_ID,z_value,z_width,N_skewers,n,k1,C0,C1,C2,beta,D0,D1,D2,pixels=[]):

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

        self.mean_F = None
        self.k_kms = None
        self.Pk_kms = None
        self.bias_delta = None
        self.bias_eta = None

        return

    #Get the values of tau0 (FGPA) for the measurement's parameter values.
    def get_tau0(self,z,z0=3.0):
        x = (1 + z)/(1 + z0)
        tau0 = np.exp(utils.quadratic_log(x,self.C0,self.C1,self.C2))
        return tau0

    #Get the values of the exponent (FGPA) for the measurement's parameter values.
    def get_exponent(self,z,z0=3.0):
        x = (1 + z)/(1 + z0)
        exponent = np.exp(utils.quadratic_log(x,self.D0,self.D1,self.D2))
        return exponent

    def get_details(self):
        details = (self.z_value,self.z_width,self.N_skewers,
                self.n,self.k1,self.C0,self.C1,self.C2,self.D0,self.D1,self.D2,self.beta,self.pixels)
        return details

    def add_Pk1D_measurement(self,pixel_object,mean_F=None):
        #Produce delta flux skewers.
        if mean_F is None:
            if not self.mean_F:
                self.add_mean_F_measurement(pixel_object)
            mean_F = self.mean_F
        F = pixel_object.lya_absorber.transmission()
        delta_F = F/mean_F - 1

        #Extract additional information from data object.
        IVAR = pixel_object.IVAR_rows
        R = pixel_object.R
        dr_hMpc = (R[-1] - R[0])/(R.shape[0] - 1)
        z = pixel_object.Z

        #Get Pk1D and add it to the measurement object.
        k_kms, Pk_kms, var_kms = Pk1D.get_Pk1D(delta_F,IVAR,dr_hMpc,z,z_value=self.z_value,z_width=self.z_width)
        self.k_kms = k_kms
        self.Pk_kms = Pk_kms

        return

    def add_mean_F_measurement(self,pixel_object):
        #Get mean_F and add it to the measurement object.
        self.mean_F = np.average(pixel_object.get_mean_quantity('flux',z_value=self.z_value,z_width=self.z_width,single_value=False))
        return

    def add_bias_delta_measurement(self,pixel_object,d=0.001,weights=None):
        #Get bias_delta and add it to the measurement object.
        self.bias_delta = bias.get_bias_delta(pixel_object,self.z_value,z_width=self.z_width,d=d,weights=weights)
        return

    def add_bias_eta_measurement(self,pixel_object,weights_dict=None,d=0.0,thermal=False,lambda_buffer=100.):
        #Get bias_eta and add it to the measurement object.
        self.bias_eta = bias.get_bias_eta(pixel_object,self.z_value,weights_dict=weights_dict,d=d,z_width=self.z_width,include_thermal_effects=thermal,lambda_buffer=lambda_buffer)
        return

    def add_sigma_dF_measurement(self,pixel_object):
        #Get sigma_dF and add it to the measurement object.
        self.sigma_dF = pixel_object.get_sigma_dF(pixel_object.lya_absorber,z_value=self.z_value,z_width=self.z_width)
        return

    def add_Pk1D_chi2(self,min_k=None,max_k=None,denom="krange",eps=0.1):
        #Calculate the "model" values.
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

        #Determine the denominator (i.e. weighting scheme).
        A = 10**10
        if denom == "uniform":
            denom = (eps * model_Pk_kms)**2
        elif denom == "krange":
            eps = A * np.ones(self.k_kms.shape)
            eps[min_j:max_j] *= eps / A
            denom = (eps * model_Pk_kms)**2
        elif denom == "krange_smooth":
            eps = A * np.ones(self.k_kms.shape)
            eps[min_j:max_j] *= eps / A
            smooth_width = (max_k - min_k)/0.01
            lower_smooth = A + (eps - A)*np.exp(-((self.k_kms - min_k)**2)/(2*smooth_width**2))
            upper_smooth = A + (eps - A)*np.exp(-((self.k_kms - max_k)**2)/(2*smooth_width**2))
            eps[:min_j] = lower_smooth[:min_j]
            eps[max_j:] = upper_smooth[max_j:]
            denom = (eps * model_Pk_kms)**2
        elif denom == "npower":
            k0 = max_k
            n = 2.
            eps = eps * ((1 + (self.k_kms/k0)**n))
            denom = (eps * model_Pk_kms)**2
        elif denom == "npower_cutoff":
            k0 = max_k
            n = 2.
            cutoff = 0.02 #kms
            eps = eps * ((1 + (self.k_kms/k0)**n))
            eps[self.k_kms>cutoff] = A
            denom = (eps * model_Pk_kms)**2

        #Add the weighting scheme and chi2 to the measurment object.
        self.Pk_kms_chi2_eps = eps
        chi2 = np.sum(((self.Pk_kms - model_Pk_kms)**2)/denom)
        self.Pk_kms_chi2 = chi2

        return

    def add_mean_F_chi2(self,eps=0.1,model='Becker13'):
        #Calculate the "model" value and weight.
        model_mean_F = get_mean_F_model(self.z_value,model=model)
        denom = (eps * model_mean_F)**2

        #Add the chi2 to the measurment object.
        self.mean_F_chi2 = np.sum(((self.mean_F - model_mean_F)**2)/denom)
        return

    def add_sigma_dF_chi2(self,min_k=None,max_k=None,eps=0.1,l_hMpc=0.25):
        #Calculate the "model" value and weight.
        model_sigma_F = get_sigma_dF_P1D(self.z_value,l_hMpc=l_hMpc)
        denom = (eps * model_sigma_F)**2

        #Add the chi2 to the measurment object.
        self.sigma_F_chi2 = np.sum(((self.sigma_F - model_sigma_F)**2)/denom)
        return

    def add_bias_delta_chi2(self,eps=0.1,model='BOSS_DR12_joint'):
        #Calculate the "model" value and weight.
        model_bias_delta, model_bias_eta = get_model_biases(self.z_value,model=model)
        denom = (eps * model_bias_delta)**2

        #Add the chi2 to the measurment object.
        self.bias_delta_chi2 = np.sum(((self.bias_delta - model_bias_delta)**2)/denom)
        return

    def add_bias_eta_chi2(self,eps=0.1,model='BOSS_DR12_joint'):
        #Calculate the "model" value and weight.
        model_bias_delta, model_bias_eta = get_model_biases(self.z_value,model=model)
        denom = (eps * model_bias_eta)**2

        #Add the chi2 to the measurment object.
        self.bias_eta_chi2 = np.sum(((self.bias_eta - model_bias_eta)**2)/denom)
        return

    def add_all_measurements_chi2(self):

        return

    @classmethod
    def combine_measurements(cls,m1,m2):
        #Check that the measurements are combinable.
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
        if utils.confirm_identical(m1.k_kms,m2.k_kms,item_name='k_kms',array=True):
            k_kms = m1.k_kms

        #Combine the different elements.
        N_skewers = m1.N_skewers + m2.N_skewers
        pixels = m1.pixels + m2.pixels

        #Combine measurements by averaging.
        mean_F = (m1.mean_F*m1.N_skewers + m2.mean_F*m2.N_skewers)/(m1.N_skewers + m2.N_skewers)
        Pk_kms = (m1.Pk_kms*m1.N_skewers + m2.Pk_kms*m2.N_skewers)/(m1.N_skewers + m2.N_skewers)
        #sigma_dF = np.sqrt(((m1.sigma_dF**2)*m1.N_skewers + (m2.sigma_dF**2)*m2.N_skewers)/(m1.N_skewers + m2.N_skewers))
        bias_delta = (m1.bias_delta*m1.N_skewers + m2.bias_delta*m2.N_skewers)/(m1.N_skewers + m2.N_skewers)
        if m1.bias_eta is not None and m2.bias_eta is not None:
            bias_eta = (m1.bias_eta*m1.N_skewers + m2.bias_eta*m2.N_skewers)/(m1.N_skewers + m2.N_skewers)
        elif m1.bias_eta is not None:
            bias_eta = m1.bias_eta
        elif m2.bias_eta is not None:
            bias_eta = m2.bias_eta
        else:
            bias_eta = None

        #Make the combined object.
        combined = cls(parameter_ID,z_value,z_width,N_skewers,n,k1,C0,C1,C2,beta,D0,D1,D2,pixels=pixels)

        #Add details to the combined object.
        combined.mean_F = mean_F
        combined.k_kms = k_kms
        combined.Pk_kms = Pk_kms
        #combined.sigma_dF = sigma_dF
        combined.bias_delta = bias_delta
        combined.bias_eta = bias_eta

        return combined

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
        dtype = [('k_kms', 'f4'), ('Pk_kms', 'f4')]
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

    def ID_filter(self,ID):
        filtered_measurements = []
        for m in self.measurements:
            if m.parameter_ID == ID:
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
        #Determine which z values we have in the set.
        z_values = list(set([m.z_value for m in self.measurements]))

        combined_measurements = []
        for z_value in z_values:
            #For each z value, filter the measurements, and get the parameter IDs.
            z_set = self.z_filter(z_value)
            z_parameter_IDs = list(set([m.parameter_ID for m in z_set.measurements]))

            for parameter_ID in z_parameter_IDs:
                #For each parameter ID, filter the measurements and combine into 1.
                z_parameter_set = z_set.ID_filter(parameter_ID)
                c_m = z_parameter_set.measurements[0]
                if len(z_parameter_set.measurements) > 1:
                    for m in z_parameter_set.measurements[1:]:
                        c_m = function_measurement.combine_measurements(c_m,m)
                combined_measurements += [c_m]

        return measurement_set(measurements=combined_measurements)

    #Don't think this is used?
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
def get_mean_F_model(z,model='Becker13'):
    if model == 'Becker13':
        mean_F = np.exp(-0.751*(((1+z)/4.5)**2.9)+0.132)
    elif model == 'FontRibera12':
        mean_F = np.exp((np.log(0.8))*(((1+z)/3.25)**3.2))
    return mean_F

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


def get_model_biases(z,model='BOSS_DR12_joint'):

    if model == 'BOSS_DR12_joint':
        #BOSS DR12 combined auto+cross data, from du Mas des Bourboux et all (2017)
        data_z = np.array([2.4])
        data_beta = np.array([1.650])
        data_beta_err = np.array([0.081])
        data_bias_delta_1plusbeta = np.array([-0.3544])
        data_bias_delta_1plusbeta_err = np.array([0.0038])

        data_f = 0.9625
        data_bias_delta_z_evol_exponent = 2.9
        # TODO: Improve on this - need to calculate f with the same parameters as in the relevant paper.
        data_bias_eta_z_evol_exponent = 2.9

        data_bias_delta = data_bias_delta_1plusbeta / (1 + data_beta)
        data_bias_eta = data_beta * data_bias_delta / data_f

    elif model == 'eBOSS_DR14_joint':
        #eBOSS DR14 combined auto+cross data, from de Saintie Agathe et al. (2019)
        data_z = np.array([2.34])
        data_beta = np.array([1.994])
        data_beta_err = np.array([0.099])
        data_bias_eta = np.array([-0.214])
        data_bias_eta_err = np.array([0.004])

        data_f = 0.96612
        data_bias_delta_z_evol_exponent = 2.9
        # TODO: Improve on this - need to calculate f with the same parameters as in the relevant paper.
        data_bias_eta_z_evol_exponent = 2.9

        data_bias_delta = (data_f * data_bias_eta) / data_beta

    z_evol_bias_delta = ((1 + z)/(1 + data_z))**data_bias_delta_z_evol_exponent
    bias_delta_model = data_bias_delta * z_evol_bias_delta

    z_evol_bias_eta = ((1 + z)/(1 + data_z))**data_bias_eta_z_evol_exponent
    bias_eta_model = data_bias_eta * z_evol_bias_eta

    return bias_delta_model, bias_eta_model


################################################################################
"""
Below: old tuning, no-RSD, theoretical based Method
"""

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
def density_to_flux(sigma_G,alpha,beta,D,mean_only=False,int_lim_fac=10.0):

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

# Not yet implemented.
"""
def load_tuning(f,mode='parameters'):

    tuning_dict = {}

    #Open the tuning file and extract the lognormal/FGPA transformation parameters.
    h = fits.open(tuning_file)
    tuning_dict['z'] = h[1].data['z']
    try:
        tuning_dict['tau0_of_z'] = h[1].data['tau0_of_z']
    except:
        tuning_dict['tau0_of_z'] = h[1].data['alpha']
    try:
        tuning_dict['texp_of_z'] = h[1].data['texp_of_z']
    except:
        tuning_dict['texp_of_z'] = h[1].data['beta']
    try:
        tuning_dict['seps_of_z'] = h[1].data['seps_of_z']
    except:
        tuning_dict['seps_of_z'] = h[1].data['sigma_G']

    #Extract additional parameters from the file's header.
    try:
        tuning_dict['tau0_A0'] = h[1].header['C0']
    except:
        tuning_dict['tau0_A0'] = h[1].header['C0']
    try:
        tuning_dict['tau0_A1'] = h[1].header['C1']
    except:
        tuning_dict['tau0_A1'] = h[1].header['C1']
    try:
        tuning_dict['tau0_A2'] = h[1].header['C2']
    except:
        tuning_dict['tau0_A2'] = h[1].header['C2']

    try:
        tuning_dict['texp_A0'] = h[1].header['D0']
    except:
        tuning_dict['texp_A0'] = h[1].header['D0']
    try:
        tuning_dict['texp_A1'] = h[1].header['D1']
    except:
        tuning_dict['texp_A1'] = h[1].header['D1']
    try:
        tuning_dict['texp_A2'] = h[1].header['D2']
    except:
        tuning_dict['texp_A2'] = h[1].header['D2']

    try:
        tuning_dict['seps_A0'] = h[1].header['D0']
    except:
        tuning_dict['seps_A0'] = h[1].header['D0']
    try:
        tuning_dict['seps_A1'] = h[1].header['D1']
    except:
        tuning_dict['seps_A1'] = h[1].header['D1']
    try:
        tuning_dict['seps_A2'] = h[1].header['D2']
    except:
        tuning_dict['seps_A2'] = h[1].header['D2']

    tuning_dict['n'] = h[1].header['n']
    tuning_dict['k1'] = h[1].header['k1']
    try:
        tuning_dict['R_kms'] = h[1].header['R_kms']
    except:
        tuning_dict['R_kms'] = h[1].header['R']
    try:
        tuning_dict['a_v'] = h[1].header['a_v']
    except:
        tuning_dict['a_v'] = h[1].header['vb']

    #Close the file.
    h.close()

    return tuning_dict

def iminuit_input_from_tuning_file(filepath):

    td = load_tuning(filepath)

    a_kwargs = {'C0' : td['tau0_A0'],      'error_C0' : 1.0,   'fix_C0' : fix_all|fix_C0,     'limit_C0' : (0., 100.),
                'C1' : td['tau0_A1'],      'error_C1' : 1.0,   'fix_C1' : fix_all|fix_C1,     #'limit_C1' : (0., 20.),
                'C2' : td['tau0_A2'],      'error_C2' : 1.0,   'fix_C2' : fix_all|fix_C2,     #'limit_C2' : (0., 20.),
                }

    b_kwargs = {'beta' : td['texp_A0'],  'error_beta' : 1.0, 'fix_beta' : fix_all|fix_beta, 'limit_beta' : (0.,5.)
                }

    sG_kwargs = {'D0' : td['seps_A0'],     'error_D0' : 1.0,   'fix_D0' : fix_all|fix_D0,     'limit_D0' : (0., 100.),
                 'D1' : td['seps_A1'],     'error_D1' : 0.2,   'fix_D1' : fix_all|fix_D1,     #'limit_D1' : (0., 20.),
                 'D2' : td['seps_A2'],     'error_D2' : 1.0,   'fix_D2' : fix_all|fix_D2,     #'limit_D2' : (0., 20.),
                 }

    s_kwargs = {'n'  : td['n'],       'error_n' : 1.0,    'fix_n' : fix_all|fix_n,       'limit_n' : (-2., 10.),
                'k1' : td['k1'],      'error_k1' : 0.001, 'fix_k1' : fix_all|fix_k1,     'limit_k1' : (0., 0.1),
                }

    other_kwargs = {'R'  : td,    'error_R' : 1.0,   'fix_R' : fix_all|fix_R,       'limit_R' : (0., 1000.),
                    'a_v': td,  'error_a_v' : 0.1, 'fix_a_v' : fix_all|fix_a_v,   'limit_a_v' : (0., 2.0),
                    'return_measurements'  : False,    'fix_return_measurements' : True,
                    'errordef'             : 1,
                    }


    return iminuit_initial
"""
