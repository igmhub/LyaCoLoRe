import numpy as np
import copy
import time

from pyacolore import utils

lya = utils.lya_rest

#Function to get the bias of delta at various z values from a sim data object.
#Assumes that the object already has tau calculated and RSDs applied.
def get_bias_delta(data,z_values,weights=None,d=0.001,z_width=0.2):

    betas = data.transformation.get_texp(data.Z)

    if isinstance(z_values, float):
        z_values = np.array([z_values])

    if not weights:
        print('no weights dict provided to get_bias_delta. Calculating weights...')
        data.compute_RSD_weights()

    #Add small extra delta to Gaussian skewers to simulate overdensity
    overdensity = copy.deepcopy(data)
    #overdensity.GAUSSIAN_DELTA_rows += d
    #overdensity.compute_physical_skewers()
    #overdensity.compute_all_tau_skewers()
    overdensity.lya_absorber.tau *= np.exp(betas*overdensity.D*d)
    overdensity.add_all_RSDs()

    #Subtract small extra delta to Gaussian skewers to simulate underdensity
    underdensity = copy.deepcopy(data)
    #underdensity.GAUSSIAN_DELTA_rows -= d
    #underdensity.compute_physical_skewers()
    #underdensity.compute_all_tau_skewers()
    underdensity.lya_absorber.tau /= np.exp(betas*underdensity.D*d)
    underdensity.add_all_RSDs()

    #Don't think this is quite valid as D(r_j) != D(s_j)
    #i.e. RSDs muddle the cells in each skewer and so this scaling would need to be done before adding RSDs
    #Same applies to the overdensity, of course
    #underdensity.lya_absorber.tau /= np.exp(betas*underdensity.D*d)

    #Calculate mean fluxes in under and overdensities, as well as normal
    biases = []
    for z_value in z_values:

        #We get means across the z-chunk and combine once bias has been computed.
        #This avoids overweighting the low end of the chunk.
        mean_F_over = overdensity.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)
        mean_F_under = underdensity.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)
        mean_F = data.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)

        #Get relevant values of D to scale by
        #(we have added d to the Gaussian field, not the density field)
        i_lower = np.searchsorted(data.Z,z_value - z_width/2.)
        i_upper = np.searchsorted(data.Z,z_value + z_width/2.)
        D_values = data.D[i_lower:i_upper]

        bias = np.average((1/mean_F) * (1/(2.*d*D_values)) * (mean_F_over - mean_F_under))
        biases += [bias]

    biases = np.array(biases)

    return biases

#Function to get the different RSD weights for calculating bias_eta.
def get_bias_eta_weights(data,z_values,d=0.0,z_width=0.2,include_thermal_effects=False,lambda_buffer=100.):

    weights_dict = {}

    for z_value in z_values:

        z_val_weights_dict = {}

        data_copy_z_val_upper = copy.deepcopy(data)
        data_copy_z_val_lower = copy.deepcopy(data)

        z_min = z_value - 0.5*z_width
        z_max = z_value + 0.5*z_width
        lambda_min = lya * (1 + z_min)
        lambda_max = lya * (1 + z_max)

        data_copy_z_val_upper.trim_skewers(lambda_min-lambda_buffer,lambda_max=lambda_max+lambda_buffer,extra_cells=1)
        data_copy_z_val_lower.trim_skewers(lambda_min-lambda_buffer,lambda_max=lambda_max+lambda_buffer,extra_cells=1)

        data_copy_z_val_upper.compute_RSD_weights(thermal=include_thermal_effects,d=d,z_r0=z_value)
        data_copy_z_val_lower.compute_RSD_weights(thermal=include_thermal_effects,d=-d,z_r0=z_value)

        z_val_weights_dict['grad_increase'] = data_copy_z_val_upper.RSD_weights
        z_val_weights_dict['grad_decrease'] = data_copy_z_val_lower.RSD_weights

        weights_dict[z_value] = z_val_weights_dict

    return weights_dict

#Function to get the bias of eta at various z values from a sim data object.
#Assumes that the object already has tau calculated but with no RSDs applied.
def get_bias_eta(data,z_values,weights_dict=None,d=0.0,z_width=0.2,include_thermal_effects=False,lambda_buffer=100.):

    alphas = data.transformation.get_tau0(data.Z)
    betas = data.transformation.get_texp(data.Z)

    if isinstance(z_values, float):
        z_values = np.array([z_values])

    """
    #Method 1:
    #Use mean of FlnF
    #data_noRSDs = copy.deepcopy(data)
    #data_noRSDs.compute_tau_skewers(data_noRSDs.lya_absorber,alphas,betas)

    biases = []
    for z_value in z_values:
        mean_FlnF = np.average(data.get_mean_quantity('FlnF',z_value=z_value,z_width=z_width,single_value=False,power=1))
        #mean_FlnFlnF = np.average(data.get_mean_quantity('FlnFlnF',z_value=z_value,z_width=z_width,single_value=False,power=1))
        mean_F = np.average(data.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1))
        biases += [np.average((mean_FlnF)/mean_F)]
    biases = np.array(biases)
    """

    """
    #Method 2:
    #Scale tau by *(1+d)

    #Add small extra grad to velocity skewers.
    grad_increase = copy.deepcopy(data)
    grad_increase.R *= (1-d)
    grad_increase.lya_absorber.tau *= (1+d)

    #Subtract small extra grad to velocity skewers.
    grad_decrease = copy.deepcopy(data)
    grad_decrease.R *= (1+d)
    grad_decrease.lya_absorber.tau *= (1-d)

    #Calculate mean fluxes in under and overdensities, as well as normal
    biases = []
    for z_value in z_values:

        #We get means across the z-chunk and combine once bias has been computed.
        #This avoids overweighting the low end of the chunk.
        mean_F_increase = grad_increase.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)
        mean_F_decrease = grad_decrease.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)
        mean_F = data.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)

        i_lower = np.searchsorted(data.Z,z_value - z_width/2.)
        i_upper = np.searchsorted(data.Z,z_value + z_width/2.)
        D_values = data.D[i_lower:i_upper]

        bias = np.average((1/mean_F) * (1/(2.*d)) * (mean_F_increase - mean_F_decrease))
        biases += [bias]

    biases = np.array(biases)
    """

    """
    #Method 3:
    #Scale tau by /(1-d)

    #Add small extra grad to velocity skewers.
    grad_increase = copy.deepcopy(data)
    grad_increase.R *= (1-d)
    grad_increase.lya_absorber.tau /= (1-d)

    #Subtract small extra grad to velocity skewers.
    grad_decrease = copy.deepcopy(data)
    grad_decrease.R *= (1+d)
    grad_decrease.lya_absorber.tau /= (1+d)

    #Calculate mean fluxes in under and overdensities, as well as normal
    biases = []
    for z_value in z_values:

        #We get means across the z-chunk and combine once bias has been computed.
        #This avoids overweighting the low end of the chunk.
        mean_F_increase = grad_increase.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)
        mean_F_decrease = grad_decrease.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)
        mean_F = data.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)

        i_lower = np.searchsorted(data.Z,z_value - z_width/2.)
        i_upper = np.searchsorted(data.Z,z_value + z_width/2.)
        D_values = data.D[i_lower:i_upper]

        bias = np.average((1/mean_F) * (1/(2.*d)) * (mean_F_increase - mean_F_decrease))
        biases += [bias]

    biases = np.array(biases)
    """

    #Method 4:
    #Whilst moving cells in RSD, add in a shift in r

    if not weights_dict:
        print('no weights dict provided to get_bias_eta. Calculating weights...')
        weights_dict = get_bias_eta_weights(data,z_values,d=d,z_width=z_width,include_thermal_effects=include_thermal_effects,lambda_buffer=lambda_buffer)

    #Copy the data and overwrite the tau skewers to remove RSDs.
    data_noRSDs = copy.deepcopy(data)
    data_noRSDs.lya_absorber.tau = data_noRSDs.lya_absorber.tau_noRSD

    #Calculate mean fluxes in under and overdensities, as well as normal
    biases = []
    for z_value in z_values:

        weights_grad_increase = weights_dict[z_value]['grad_increase']
        weights_grad_decrease = weights_dict[z_value]['grad_decrease']

        z_min = z_value - 0.5*z_width
        z_max = z_value + 0.5*z_width
        lambda_min = lya * (1 + z_min)
        lambda_max = lya * (1 + z_max)

        #Copy the data and then trim it to the area around the z value.
        data_noRSDs_z_val = copy.deepcopy(data_noRSDs)
        data_noRSDs_z_val.trim_skewers(lambda_min-lambda_buffer,lambda_max=lambda_max+lambda_buffer,extra_cells=1)

        #Copy the noRSD data, add RSDs with the extra shift (increase), then trim the skewers.
        grad_increase = copy.deepcopy(data_noRSDs_z_val)
        grad_increase.add_all_RSDs(weights=weights_grad_increase,thermal=include_thermal_effects)
        grad_increase.trim_skewers(lambda_min,lambda_max=lambda_max)

        #Copy the noRSD data, add RSDs with the extra shift (decrease), then trim the skewers.
        grad_decrease = copy.deepcopy(data_noRSDs_z_val)
        grad_decrease.add_all_RSDs(weights=weights_grad_decrease,thermal=include_thermal_effects)
        grad_decrease.trim_skewers(lambda_min,lambda_max=lambda_max)

        #Make small pixel object with RSDs and trim it. Also trim the noRSD object.
        data_z_val = copy.deepcopy(data)
        data_z_val.trim_skewers(lambda_min,lambda_max=lambda_max)
        data_noRSDs_z_val.trim_skewers(lambda_min,lambda_max=lambda_max)

        #We get means across the z-chunk and combine once bias has been computed.
        #This avoids overweighting the low end of the chunk.
        mean_F_increase = grad_increase.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)
        mean_F_decrease = grad_decrease.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)
        mean_F = data_z_val.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)

        bias = np.average((1/mean_F) * (1/(2.*d)) * (mean_F_increase - mean_F_decrease))
        biases += [bias]

    biases = np.array(biases)

    """
    z_min = np.min(z_values) - 0.5*z_width
    z_max = np.max(z_values) + 0.5*z_width
    lambda_min = lya * (1 + z_min)
    lambda_max = lya * (1 + z_max)

    lambda_buffer = 100. #A
    min_catalog_z = 1.8

    #Copy the data and then trim it to the area around the z value.
    data_noRSDs = copy.deepcopy(data)
    data_noRSDs.compute_tau_skewers(data_noRSDs.lya_absorber)

    #Copy the noRSD data, add RSDs with the extra shift (increase), then trim the skewers.
    grad_increase = copy.deepcopy(data_noRSDs)
    grad_increase.add_all_RSDs(thermal=include_thermal_effects,d=d,z_r0=z_r0)
    grad_increase.trim_skewers(lambda_min,min_catalog_z,lambda_max=lambda_max)

    #Copy the noRSD data, add RSDs with the extra shift (decrease), then trim the skewers.
    grad_decrease = copy.deepcopy(data_noRSDs)
    grad_decrease.add_all_RSDs(thermal=include_thermal_effects,d=-d,z_r0=z_r0)
    grad_decrease.trim_skewers(lambda_min,min_catalog_z,lambda_max=lambda_max)
    biases = []

    for z_value in z_values:
        #We get means across the z-chunk and combine once bias has been computed.
        #This avoids overweighting the low end of the chunk.
        mean_F_increase = grad_increase.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)
        mean_F_decrease = grad_decrease.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)
        mean_F = data.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)

        bias = np.average((1/mean_F) * (1./(2.*d)) * (mean_F_increase - mean_F_decrease))
        biases += [bias]
    biases = np.array(biases)
    """

    return biases


#Function to get beta at various z values from a sim data object.
#Assumes that the object already has tau calculated but with no RSDs applied.
def get_beta(data,alphas,betas,f,z_values,d=0.001,z_width=0.2,z_r0=2.5,include_thermal_effects=False):

    biases_delta = get_bias_delta(data,betas,z_values,d=d,z_width=z_width)
    biases_nu = get_bias_nu(data,alphas,betas,z_values,d=d,z_width=z_width,z_r0=z_r0,include_thermal_effects=include_thermal_effects)

    betas_RSD = f*biases_nu/biases_delta

    return betas_RSD
