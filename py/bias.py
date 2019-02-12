import numpy as np
import copy

#Function to get the bias of delta at various z values from a sim data object.
#Assumes that the object already has tau calculated and RSDs applied.
def get_bias_delta(data,z_values,d=0.05,z_width=0.2,RSD_weights=None):

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

    #Calculate mean fluxes in under and overdensities, as well as normal
    mean_F_over = []
    mean_F_under = []
    mean_F = []
    for z_value in z_values:
        mean_F_over += [overdensity.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=True,power=1)]
        mean_F_under += [underdensity.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=True,power=1)]
        mean_F += [data.get_mean_quantity('flux',z_value=z_value,z_width=z_width,single_value=True,power=1)]

    mean_F_over = np.array(mean_F_over)
    mean_F_under = np.array(mean_F_under)
    mean_F = np.array(mean_F)

    #Calculate bias using gradient
    b = (1/mean_F) * (1/d) * (mean_F_over - mean_F_under

    return b
