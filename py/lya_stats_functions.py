import numpy as np

#Code to analyse the statistics of Lya forest data

#Produce a binned average density skewer
def density_stats(z,z_qso,density_skewers):
    """Determine the average density, density variance and delta variance over a number of redshift bins"""

    N_cells = density_skewers.shape[1]
    density_squared_skewers = density_skewers**2
    N_qso = len(z_qso)

    #Make a mask to remove irrelevant data when averaging (densities for z>z_qso)
    mask = np.ones(density_skewers.shape)
    max_pixel_qso = [0.]*N_qso
    for j in range(N_qso):
        max_pixel_qso[j] = (np.argmax(z>z_qso[j]))%N_cells
        mask[j,max_pixel_qso[j]+1:]=np.zeros(1,(mask.shape[1]-max_pixel_qso[j]-1))

    #Calculate the average density over all skewers
    mean_density = np.zeros(N_cells)
    ignore = []
    for j in range(N_cells):
        if sum(mask[:,j])!=0:
            mean_density[j] = np.average(density_skewers[:,j],weights=mask[:,j])
        else:
            ignore += [j]

    #Calculate the average density squared over all skewers
    mean_density_squared = np.zeros(N_cells)
    for j in range(N_cells):
        if sum(mask[:,j])!=0:
            mean_density_squared[j] = np.average(density_squared_skewers[:,j],weights=mask[:,j])

    #Set up the binned data structure
    N_bins = 100
    binned_z = np.zeros(N_bins)
    binned_z_location = np.zeros(N_bins)
    binned_mean_density = np.zeros(N_bins)
    binned_mean_density_alt = np.zeros(N_bins)
    binned_density_var = np.zeros(N_bins)
    binned_delta_var = np.zeros(N_bins)

    if N_cells<N_bins:
        #ERROR, not enough cells for the N_bins specified
        exit('Not enough cells for the N_bins specified.')

    #Construct the conversion matrix: an N_bins by N_cells matrix detailing the contribution of each cell to each bin
    #conversion[i,j] denotes the contribution of the jth cell to the ith bin
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
        #Construct output vectors
        binned_z_location[i] = (i+0.5)*ratio + 0.5

        binned_mean_density[i] = sum(conversion[i,:]*mean_density)/ratio
        binned_density_var[i] = sum(conversion[i,:]*mean_density_squared)/ratio - (binned_mean_density[i])**2
        binned_delta_var[i] = binned_density_var[i]/((binned_mean_density[i])**2)

    #Calculate the binned z values by linear interpolation
    binned_z = np.interp(binned_z_location,list(range(N_cells)),z)

    #for j in range(N_bins):
    #    conversion_mask = mask*conversion[j,:]
    #    binned_mean_density_alt[j] = sum(sum(density_skewers*conversion_mask))/sum(sum(conversion_mask))

    return binned_z, binned_mean_density, binned_density_var, binned_delta_var
