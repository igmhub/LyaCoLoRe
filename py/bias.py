import numpy as np
import copy

#Function to get the bias of delta at various z values from a sim data object.
#Assumes that the object already has tau calculated and RSDs applied.
def get_bias_delta(data,betas,z_values,d=0.001,z_width=0.2):

    """
    #Add small extra delta to Gaussian skewers to simulate overdensity
    overdensity = copy.deepcopy(data)
    overdensity.GAUSSIAN_DELTA_rows += d
    overdensity.compute_physical_skewers()
    overdensity.compute_all_tau_skewers(alphas,betas)
    overdensity.add_all_RSDs(thermal=False,weights=RSD_weights)

    #Subtract small extra delta to Gaussian skewers to simulate underdensity
    underdensity = copy.deepcopy(data)
    underdensity.GAUSSIAN_DELTA_rows -= d
    underdensity.compute_physical_skewers()
    underdensity.compute_all_tau_skewers(alphas,betas)
    underdensity.add_all_RSDs(thermal=False,weights=RSD_weights)
    """

    #Add small extra delta to Gaussian skewers to simulate overdensity
    overdensity = copy.deepcopy(data)
    overdensity.GAUSSIAN_DELTA_rows += d
    overdensity.DENSITY_DELTA_rows = (overdensity.DENSITY_DELTA_rows + 1)*np.exp(overdensity.D*d) - 1
    overdensity.lya_absorber.tau *= np.exp(betas*overdensity.D*d)

    #Subtract small extra delta to Gaussian skewers to simulate underdensity
    underdensity = copy.deepcopy(data)
    underdensity.GAUSSIAN_DELTA_rows -= d
    underdensity.DENSITY_DELTA_rows = (underdensity.DENSITY_DELTA_rows + 1)*np.exp(underdensity.D*d) - 1
    underdensity.lya_absorber.tau /= np.exp(betas*underdensity.D*d)

    #Calculate mean fluxes in under and overdensities, as well as normal
    biases = []
    for z_value in z_values:

        #We get means across the z-chunk and combine once bias has been computed.
        #This avoids overweighting the low end of the chunk.
        mean_F_over = overdensity.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)
        mean_F_under = underdensity.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)
        mean_F = data.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=False,power=1)

        i_lower = np.searchsorted(data.Z,z_value - z_width/2.)
        i_upper = np.searchsorted(data.Z,z_value + z_width/2.)
        D_values = data.D[i_lower:i_upper]
        #print(np.sum(data.IVAR_rows[:,i_lower:i_lower+5],axis=0),np.sum(data.IVAR_rows[:,i_upper-5:i_upper],axis=0))

        #D_value = np.interp(z_value,data.Z,data.D)
        bias = np.average((1/mean_F) * (1/(2.*d*D_values)) * (mean_F_over - mean_F_under))
        biases += [bias]

    biases = np.array(biases)

    return biases

#Function to get the different RSD weights for calculating bias_nu.


#Function to get the bias of eta at various z values from a sim data object.
#Assumes that the object already has tau calculated but with no RSDs applied.
def get_bias_nu(data,alphas,betas,z_values,d=0.001,z_width=0.2,z_r0=2.5,include_thermal_effects=False):

    #print('orig RSD:\n',data.lya_absorber.tau[0,:5])
    #Copy the data and overwrite the tau skewers to remove RSDs.
    data_noRSDs = copy.deepcopy(data)
    data_noRSDs.compute_tau_skewers(data_noRSDs.lya_absorber,alphas,betas)
    #print('noRSD:\n',data_noRSDs.lya_absorber.tau[0,:5])

    #Add small extra grad to velocity skewers.
    grad_increase = copy.deepcopy(data_noRSDs)
    grad_increase.add_all_RSDs(thermal=include_thermal_effects,d=d,z_r0=z_r0)
    #print('new RSD:\n',grad_increase.lya_absorber.tau[0,:5])

    #Subtract small extra grad to velocity skewers.
    grad_decrease = copy.deepcopy(data_noRSDs)
    grad_decrease.add_all_RSDs(thermal=include_thermal_effects,d=-d,z_r0=z_r0)
    #print('new RSD:\n',grad_decrease.lya_absorber.tau[0,:5])

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
        #print(np.sum(data.IVAR_rows[:,i_lower:i_lower+5],axis=0),np.sum(data.IVAR_rows[:,i_upper-5:i_upper],axis=0))

        #D_value = np.interp(z_value,data.Z,data.D)
        #bias = np.average((1/mean_F) * (1/(2.*d*D_values)) * (mean_F_increase - mean_F_decrease))
        bias = np.average((1/mean_F) * (1/(2.*d)) * (mean_F_increase - mean_F_decrease))
        biases += [bias]

    biases = np.array(biases)

    """
    biases = []
    for z_value in z_values:
        mean_FlnF = np.average(data.get_mean_quantity('FlnF',z_value=z_value,z_width=z_width,single_value=False,power=1))
        biases += [mean_FlnF]
    biases = np.array(biases)
    """

    """
    #Add small extra grad to velocity skewers.
    grad_increase = copy.deepcopy(data)
    grad_increase.R *= (1-d)
    grad_increase.Z = np.interp(grad_increase.R,data.R,data.Z)
    grad_increase.Z_QSO = np.interp(data.Z_QSO,data.Z,grad_increase.Z)

    #Subtract small extra grad to velocity skewers.
    grad_decrease = copy.deepcopy(data)
    grad_decrease.R *= (1+d)
    grad_decrease.Z = np.interp(grad_decrease.R,data.R,data.Z)
    grad_decrease.Z_QSO = np.interp(data.Z_QSO,data.Z,grad_decrease.Z)

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
        #print(np.sum(data.IVAR_rows[:,i_lower:i_lower+5],axis=0),np.sum(data.IVAR_rows[:,i_upper-5:i_upper],axis=0))

        #D_value = np.interp(z_value,data.Z,data.D)
        bias = np.average((1/mean_F) * (1/(2.*d*D_values)) * (mean_F_increase - mean_F_decrease))
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