################################################################################

"""
We would like to add small scale flucatuations to the Gaussian field.
We must work out how much variance to add.
This is done by
 - computing the analytical P1D flux variance from Palanque-Delabrouille et al. (2013),
 - computing the Gaussian variance required to achieve this flux variance using Font-Ribera et al. (2012), and thus the extra variance our Gaussian skewers require
 - stretching the current skewers to achieve smaller cell sizes (NGP)
 - adding a random number to each cell with the appropriate statistics
"""

def tune():

    #Work out sigma_G desired to achive the P1D sigma_dF
    tuning_z_values = np.linspace(0,4.0,128)
    beta = 1.65

    print('\nCalculating how much extra power to add...')

    D_values = np.interp(tuning_z_values,cosmology_data['Z'],cosmology_data['D'])
    sigma_G_tolerance = 0.0001

    def tune_sigma_G(z,D,l_hMpc,beta,Om):

        sigma_dF_needed = functions.get_sigma_dF_P1D(z,l_hMpc=l_hMpc,Om=Om)
        mean_F_needed = functions.get_mean_F_model(z)

        alpha,sigma_G,mean_F,sigma_dF = functions.find_sigma_G(mean_F_needed,sigma_dF_needed,beta,D,tolerance=sigma_G_tolerance)

        return (z,alpha,sigma_G,mean_F,sigma_dF,mean_F_needed,sigma_dF_needed)

    tasks = [(z,np.interp(z,cosmology_data['Z'],cosmology_data['D']),final_cell_size,beta,simulation_parameters['omega_M']) for z in tuning_z_values]

    if __name__ == '__main__':
        pool = Pool(processes = N_processes)
        results = []
        start_time = time.time()

        for task in tasks:
            pool.apply_async(tune_sigma_G,task,callback=log_result,error_callback=log_error)

        pool.close()
        pool.join()

    dtype = [('z', 'f4'), ('alpha', 'f4'), ('sigma_G', 'f4'), ('mean_F', 'f4'), ('sigma_dF', 'f4'), ('mean_F_needed', 'f4'), ('sigma_dF_needed', 'f4')]
    tune_small_scale_fluctuations = np.array(results,dtype=dtype)
    tune_small_scale_fluctuations = np.sort(tune_small_scale_fluctuations,order=['z'])

    plt.plot(tune_small_scale_fluctuations['z'],tune_small_scale_fluctuations['alpha'],label='alpha')
    plt.plot(tune_small_scale_fluctuations['z'],tune_small_scale_fluctuations['sigma_G'],label='sigma_G')
    plt.plot(tune_small_scale_fluctuations['z'],tune_small_scale_fluctuations['mean_F'],label='mean_F')
    plt.plot(tune_small_scale_fluctuations['z'],tune_small_scale_fluctuations['sigma_dF'],label='sigma_dF')
    plt.grid()
    plt.legend()
    plt.savefig('tune_flux_values_tol{}_n{}.pdf'.format(sigma_G_tolerance,tuning_z_values.shape[0]))
    plt.show()
    """
    plt.plot(tune_small_scale_fluctuations['z'],tune_small_scale_fluctuations['mean_F']/tune_small_scale_fluctuations['mean_F_needed'] - 1,label='mean_F error')
    plt.plot(tune_small_scale_fluctuations['z'],tune_small_scale_fluctuations['sigma_dF']/tune_small_scale_fluctuations['sigma_dF_needed'] - 1,label='sigma_dF error')
    plt.grid()
    plt.legend()
    plt.savefig('tune_flux_values_tol{}_n{}_Ferrors.pdf'.format(sigma_G_tolerance,tuning_z_values.shape[0]))
    #plt.show()
    """
    header = fits.Header()
    header['beta'] = beta
    header['l_hMpc'] = final_cell_size

    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)
    cols_DATA = fits.ColDefs(tune_small_scale_fluctuations)
    hdu_DATA = fits.BinTableHDU.from_columns(cols_DATA,header=header,name='DATA')

    hdulist = fits.HDUList([prihdu, hdu_DATA])
    hdulist.writeto('input_files/tune_small_scale_fluctuations_n{}.fits'.format(tuning_z_values.shape[0]))
    hdulist.close()

    return tuning_z_values
