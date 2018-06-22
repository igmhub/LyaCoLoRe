import numpy as np
from astropy.io import fits

#Function to calculate the mean of deltas, mean of deltas^2, and N.
def return_means(DELTA_rows,weights,sample_pc=1.0):
    DELTA_SQUARED_rows = DELTA_rows**2
    N_cells = DELTA_rows.shape[1]

    N = np.zeros(N_cells)
    mean_DELTA = np.zeros(N_cells)
    mean_DELTA_SQUARED = np.zeros(N_cells)

    for j in range(N_cells):
        N[j] = np.sum(weights[:,j],axis=0)
        if N[j] > 0:
            mean_DELTA[j] = np.average(DELTA_rows[:,j],weights=weights[:,j])
            mean_DELTA_SQUARED[j] = np.average(DELTA_SQUARED_rows[:,j],weights=weights[:,j])

    return N, mean_DELTA, mean_DELTA_SQUARED

#
def combine_pixel_means(results):

    N_cells = results[0][0].shape[0]
    N = np.zeros(N_cells)
    mean_DELTA = np.zeros(N_cells)
    mean_DELTA_SQUARED = np.zeros(N_cells)

    for result in results:
        N += result[0]

    for j in range(N_cells):
        if N[j] > 0:
            for result in results:
                mean_DELTA[j] += result[0][j]*result[1][j]/N[j]
                mean_DELTA_SQUARED[j] += result[0][j]*result[2][j]/N[j]

    var_DELTA = mean_DELTA_SQUARED - mean_DELTA**2

    return N, mean_DELTA, var_DELTA

#Function to take a list of sets of statistics (as produced by 'get_statistics'), and calculate means and variances.
def combine_means(means_list):

    means_shape = means_list[0].shape
    means_data_type = means_list[0].dtype
    N_cells = means_shape[0]

    quantities = ['GAUSSIAN_DELTA','GAUSSIAN_DELTA_SQUARED','DENSITY_DELTA','DENSITY_DELTA_SQUARED','F','F_SQUARED','F_DELTA','F_DELTA_SQUARED']

    combined_means = np.zeros(means_shape,dtype=means_data_type)
    for means_array in means_list:
        combined_means['N'] += means_array['N']

    for i in range(N_cells):
        if combined_means['N'][i] > 0:
            for quantity in quantities:
                for means_array in means_list:
                    combined_means[quantity][i] += (means_array['N'][i]*means_array[quantity][i])/combined_means['N'][i]

    return combined_means

#Function to convert a set of means of quantities and quantities squared (as outputted by 'combine_means') to a set of means and variances.
def means_to_statistics(means):

    statistics_dtype = [('N', 'f4')
        , ('GAUSSIAN_DELTA_MEAN', 'f4'), ('GAUSSIAN_DELTA_VAR', 'f4')
        , ('DENSITY_DELTA_MEAN', 'f4'), ('DENSITY_DELTA_VAR', 'f4')
        , ('F_MEAN', 'f4'), ('F_VAR', 'f4')
        , ('F_DELTA_MEAN', 'f4'), ('F_DELTA_VAR', 'f4')]

    statistics = np.zeros(means.shape,dtype=statistics_dtype)

    statistics['N'] = means['N']
    statistics['GAUSSIAN_DELTA_MEAN'] = means['GAUSSIAN_DELTA']
    statistics['DENSITY_DELTA_MEAN'] = means['DENSITY_DELTA']
    statistics['F_MEAN'] = means['F']
    statistics['F_DELTA_MEAN'] = means['F_DELTA']

    statistics['GAUSSIAN_DELTA_VAR'] = means['GAUSSIAN_DELTA_SQUARED'] - means['GAUSSIAN_DELTA']**2
    statistics['DENSITY_DELTA_VAR'] = means['DENSITY_DELTA_SQUARED'] - means['DENSITY_DELTA']**2
    statistics['F_VAR'] = means['F_SQUARED'] - means['F']**2
    statistics['F_DELTA_VAR'] = means['F_DELTA_SQUARED'] - means['F_DELTA']**2

    return statistics

#Function to write the statistics data to file, along with an HDU extension contanint cosmology data.
def write_statistics(location,N_side,statistics,cosmology_data):

    #Construct HDU from the statistics array.
    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)
    cols_stats = fits.ColDefs(statistics)
    hdu_stats = fits.BinTableHDU.from_columns(cols_stats,name='STATISTICS')
    cols_cosmology = fits.ColDefs(cosmology_data)
    hdu_cosmology = fits.BinTableHDU.from_columns(cols_cosmology,name='COSMO')

    #Put the HDU into an HDUlist and save as a new file. Close the HDUlist.
    filename = '/statistics.fits'
    hdulist = fits.HDUList([prihdu,hdu_stats,hdu_cosmology])
    hdulist.writeto(location+filename)
    hdulist.close

    return

#Function to calculate mean_F and sigma_dF for given values of sigma_G, alpha and beta.
def get_flux_stats(sigma_G,alpha,beta,D,mean_only=False,int_lim_fac=10.0):

    int_lim = sigma_G*int_lim_fac

    delta_G_integral = np.linspace(-int_lim,int_lim,10**4)
    delta_G_integral = np.reshape(delta_G_integral,(1,delta_G_integral.shape[0]))

    prob_delta_G = (1/((np.sqrt(2*np.pi))*sigma_G))*np.exp(-(delta_G_integral**2)/(2*(sigma_G**2)))

    density_integral = gaussian_to_lognormal_delta(delta_G_integral,sigma_G,D) + 1
    F_integral = density_to_flux(density_integral,alpha,beta)

    mean_F = np.trapz(prob_delta_G*F_integral,delta_G_integral)[0]

    if mean_only == False:
        delta_F_integral = F_integral/mean_F - 1
        integrand = prob_delta_G*(delta_F_integral**2)
        sigma_dF = (np.sqrt(np.trapz(integrand,delta_G_integral)[0]))
    else:
        sigma_dF = None

    return mean_F, sigma_dF
