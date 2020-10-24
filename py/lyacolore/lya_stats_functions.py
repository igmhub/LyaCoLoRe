import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic

#Code to analyse the statistics of Lya forest data

#Function to normalise the value of delta from CoLoRe output delta skewers such that its mean is 0.
def normalise_delta(mask,delta_skewers):

    N_cells = delta_skewers.shape[1]
    mean = np.average(delta_skewers,weights=mask)
    normalised_delta_skewers = (delta_skewers - mean)/(mean + 1)

    return normalised_delta_skewers


#Function to generate a mask for an array of density skewers.
#For each qso, this denotes cells in a matrix of skewers as 'to be ignored' if z_cell > z_qso
def make_mask(z,z_qso):

    N_cells = len(z)
    N_qso = len(z_qso)

    mask = np.ones((N_qso,N_cells))
    max_pixel_qso = [0.]*N_qso

    for i in range(N_qso):
        max_pixel_qso[i] = (np.argmax(z>z_qso[i])-1)%N_cells
        mask[i,max_pixel_qso[i]+1:]=np.zeros(1,(mask.shape[1]-max_pixel_qso[i]-1))

    return mask


#Construct the conversion matrix: an N_bins by N_cells matrix detailing the contribution of each cell to each bin
#conversion[i,j] denotes the contribution of the jth cell to the ith bin
def convert_to_bins(N_bins,N_cells):

    conversion = np.zeros((N_bins,N_cells))
    ratio = N_cells/N_bins

    for i in range(N_bins):
        unassigned_row_contribution=ratio

        for j in range(N_cells):
            unassigned_col_contribution=1-sum(conversion[:,j])

            if unassigned_row_contribution>0:
                conversion[i,j] = min(unassigned_row_contribution,unassigned_col_contribution,1)
                unassigned_row_contribution -= conversion[i,j]
                unassigned_col_contribution -= conversion[i,j]

            else:
                break

    conversion = conversion/ratio

    return conversion


#Produce binned skewers of average density, density variance and delta variance.
def density_stats(z,z_qso,density_skewers):

    """Determine the average density, density variance and delta variance over a number of redshift bins"""

    #Define uesful numbers and calculate the square of density in each cell.
    N_cells = density_skewers.shape[1]
    N_qso = len(z_qso)
    density_squared_skewers = density_skewers**2

    #Define the desired number of bins.
    N_bins = 100

    #Make a mask to remove irrelevant data when averaging (densities for z>z_qso)
    mask = make_mask(z,z_qso)

    #Ensure that if we have a zero-sum column(s) in mask, it is excluded.
    #Such columns are trimmed from all matrices.
    sum_of_weights = np.sum(mask,axis=0)
    zero_weight_start = np.argmax(sum_of_weights==0)

    if zero_weight_start==0 and sum_of_weights[0]==0:
        exit()

    elif zero_weight_start>0:
        z = z[0:zero_weight_start]
        density_skewers = density_skewers[:,0:zero_weight_start]
        density_squared_skewers = density_squared_skewers[:,0:zero_weight_start]
        mask = mask[:,0:zero_weight_start]
        N_cells = zero_weight_start

    #Convert from pixel redshifts to redshift bins
    conversion = convert_to_bins(N_bins,N_cells)

    #Mask the density skewers.
    masked_density_skewers = mask*density_skewers
    masked_density_squared_skewers = mask*density_squared_skewers

    #Calculate an average density for each cell, masking irrelevant cells.
    mean_density = np.average(density_skewers,axis=0,weights=mask)
    mean_density_squared = np.average(density_squared_skewers,axis=0,weights=mask)

    #Re-bin the mean_density vector.
    binned_mean_density = np.dot(mean_density,np.transpose(conversion))
    binned_mean_density_squared = np.dot(mean_density_squared,np.transpose(conversion))

    #Calculate the variance in density and delta.
    binned_density_var = binned_mean_density_squared - binned_mean_density**2
    binned_delta_var = binned_density_var/(binned_mean_density**2)

    #Calculate the binned z values by linear interpolation
    ratio = N_cells/N_bins
    binned_z_location = np.asarray([(i+0.5)*ratio + 0.5 for i in range(N_bins)])
    binned_z = np.interp(binned_z_location,list(range(N_cells)),z)

    return binned_z, binned_mean_density, binned_density_var, binned_delta_var
