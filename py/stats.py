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
def combine_statistics(statistics_list):

    statistics_shape = statistics_list[0].shape
    statistics_data_type = statistics_list[0].dtype
    N_cells = statistics_shape[0]

    quantities = [('N', 'f4')
        , ('GAUSSIAN_DELTA_MEAN', 'f4'), ('GAUSSIAN_DELTA_VAR', 'f4')
        , ('DENSITY_DELTA_MEAN', 'f4'), ('DENSITY_DELTA_VAR', 'f4')
        , ('TAU_MEAN', 'f4'), ('TAU_VAR', 'f4')
        , ('F_MEAN', 'f4'), ('F_VAR', 'f4')]

    combined_statistics = np.zeros(statistics_shape,dtype=means_data_type)
    for s_array in statistics_list:
        combined_means['N'] += s_array['N']

    #Combine the means.
    cells = combined_means['N']>0
    for q in quantities:
        for s_array in statistics_list:
            combined_statistics[quantity][cells] += s_array['N'][cells]*s_array[q][cells]/combined_statistics['N'][cells]

    return combined_statistics

#Function to convert a set of means of quantities and quantities squared (as outputted by 'combine_means') to a set of means and variances.
def means_to_statistics(means):

    statistics_dtype = [('N', 'f4')
        , ('GAUSSIAN_DELTA_MEAN', 'f4'), ('GAUSSIAN_DELTA_VAR', 'f4')
        , ('DENSITY_DELTA_MEAN', 'f4'), ('DENSITY_DELTA_VAR', 'f4')
        , ('TAU_MEAN', 'f4'), ('TAU_VAR', 'f4')
        , ('F_MEAN', 'f4'), ('F_VAR', 'f4')]
        #, ('F_DELTA_MEAN', 'f4'), ('F_DELTA_VAR', 'f4')]

    statistics = np.zeros(means.shape,dtype=statistics_dtype)

    statistics['N'] = means['N']
    statistics['GAUSSIAN_DELTA_MEAN'] = means['GAUSSIAN_DELTA']
    statistics['DENSITY_DELTA_MEAN'] = means['DENSITY_DELTA']
    statistics['TAU_MEAN'] = means['TAU']
    statistics['F_MEAN'] = means['F']
    #statistics['F_DELTA_MEAN'] = means['F_DELTA']

    statistics['GAUSSIAN_DELTA_VAR'] = means['GAUSSIAN_DELTA_SQUARED'] - means['GAUSSIAN_DELTA']**2
    statistics['DENSITY_DELTA_VAR'] = means['DENSITY_DELTA_SQUARED'] - means['DENSITY_DELTA']**2
    statistics['TAU_VAR'] = means['TAU_SQUARED'] - means['TAU']**2
    statistics['F_VAR'] = means['F_SQUARED'] - means['F']**2
    #statistics['F_DELTA_VAR'] = means['F_DELTA_SQUARED'] - means['F_DELTA']**2

    return statistics

#Function to write the statistics data to file, along with an HDU extension contanint cosmology data.
def write_statistics(location,filename,statistics):

    #Construct HDU from the statistics array.
    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)
    cols_stats = fits.ColDefs(statistics)
    hdu_stats = fits.BinTableHDU.from_columns(cols_stats,name='STATISTICS')
    #cols_cosmology = fits.ColDefs(cosmology_data)
    #hdu_cosmology = fits.BinTableHDU.from_columns(cols_cosmology,name='COSMO')

    #Put the HDU into an HDUlist and save as a new file. Close the HDUlist.
    hdulist = fits.HDUList([prihdu,hdu_stats])
    hdulist.writeto(location+filename)
    hdulist.close

    return
