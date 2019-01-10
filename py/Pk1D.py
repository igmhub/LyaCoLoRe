import numpy as np

from . import utils

def get_Pk1D(skewer_rows,IVAR_rows,dr_hMpc,z,z_value=0.0,z_width=None,R1=25.0,units='km/s',gaussian=False):

    if z_width:
        #Find relevant chunk of the skewers
        z_min = z_value - z_width/2.
        z_max = z_value + z_width/2.

        j_lower = np.searchsorted(z,z_min)
        j_upper = np.searchsorted(z,z_max) - 1
        N_cells_chunk = j_upper - j_lower + 1
    else:
        j_lower = 0
        j_upper = skewer_rows.shape[1] - 1
        N_cells_chunk = skewer_rows.shape[1]

    #if skewer contains entire chunk, keep, otherwise discard
    N_qso = skewer_rows.shape[0]
    skewer_rows_chunk = []
    IVAR_rows_chunk = []
    for i in range(N_qso):
        if np.sum(IVAR_rows[i][j_lower:j_upper+1]) == N_cells_chunk:
            skewer_rows_chunk.append(skewer_rows[i][j_lower:j_upper+1])
            IVAR_rows_chunk.append(IVAR_rows[i][j_lower:j_upper+1])
    skewer_rows_chunk = np.array(skewer_rows_chunk)
    IVAR_rows_chunk = np.array(IVAR_rows_chunk)

    relevant_QSOs = np.sum(IVAR_rows[:,j_lower:j_upper+1],axis=1) == N_cells_chunk
    #not_relevant_QSOs = np.sum(IVAR_rows[:,j_lower:j_upper+1],axis=1) < N_cells_chunk

    #trim z to the chunk now being considered
    z = z[j_lower:j_upper+1]

    if units == 'km/s':
        #convert to kms
        #If we're dealing with the Gaussian field, we want z=0 value.
        if gaussian:
            dkms_dhMpc = utils.get_dkms_dhMpc(0.0)
        else:
            dkms_dhMpc = utils.get_dkms_dhMpc(z_value)

        #get the cell width (this is not constant in kms)
        dv_kms = dkms_dhMpc*dr_hMpc

        #get the k frequencies
        k_kms = np.fft.rfftfreq(N_cells_chunk)*2*np.pi/dv_kms

        #ft the skewers
        ft_rows = np.fft.rfft(skewer_rows_chunk,axis=1) / np.sqrt(N_cells_chunk/dv_kms)
        pk_rows = np.abs(ft_rows)**2

        #compute Fourier transform of Top-Hat filter of size l_kms and apply it
        W_kms = (np.sinc((k_kms*dv_kms)/(2*np.pi)))**2 #* np.exp(-pow(k_kms*R1,2))
        pk_rows /= (W_kms)
        #print(W_kms[-5:]**2)

        #calculate mean and variance
        pk_kms = np.average(pk_rows,axis=0)
        var_kms = np.average((pk_rows-pk_kms)**2,axis=0)

        #Relabel for return.
        pk = pk_kms
        k = k_kms
        var = var_kms

    elif units == 'Mpc/h':
        #get the k frequencies
        k_hMpc = np.fft.rfftfreq(N_cells_chunk)*2*np.pi/dr_hMpc

        #ft the skewers
        ft_rows = np.fft.rfft(skewer_rows_chunk,axis=1) / np.sqrt(N_cells_chunk/dr_hMpc)
        pk_rows = np.abs(ft_rows)**2

        #compute Fourier transform of Top-Hat filter of size l_kms and apply it
        W_hMpc = np.sinc((k_hMpc*dr_hMpc)/(2*np.pi))
        pk_rows /= (W_hMpc**2)

        #calculate mean and variance
        pk_hMpc = np.average(pk_rows,axis=0)
        var_hMpc = np.sum((pk_rows-pk_hMpc)**2,axis=0)

        #Relabel for return.
        pk = pk_hMpc
        k = k_hMpc
        var = var_hMpc

    else:
        print('Units not recognised. Please choose from \'km/s\' and \'Mpc/h\'.')

    return k, pk, var


"""
Not needed at the moment

def get_binned_separations(R_hMpc,bins,rmin=0.0,rmax=160.,nr=40):

    N_cells = R_hMpc.shape[0]

    #Calculate separations
    separations = np.zeros((N_cells,N_cells))
    for i in range(N_cells):
        for j in range(N_cells):
            separation = abs(R_hMpc[i]-R_hMpc[j])
            if separation < rmax and separation >= rmin:
                separations[i,j] = separation
            elif separation >= rmax or separation < rmin:
                separations[i,j] = -1.0

    #Put the seaparations into bins.
    binned_separations = np.digitize(separations,bins)

    #Subtract 1 from the binned_separations to start the bin numbering from 0.
    #Also subtract an identity matrix as we do not want to include a pixel's correlation with itself.
    binned_separations = binned_separations - 1 - np.identity(N_cells)

    return binned_separations


def get_cf1D(skewer_rows,R_hMpc,N_processes=1):

    rmax = 160.0
    rmin = 0.0
    nr = 40
    bins = np.linspace(rmin,rmax,nr+1)

    N_bin_splits = 100
    N_skewers = skewer_rows.shape[0]

    binned_separations = get_binned_separations(R_hMpc,bins)

    R_binned = np.zeros(nr)
    for i in range(nr):
        R_binned[i] = bins[i]+(bins[i+1]-bins[i])/2

    #For each bin, make a list of the pixel pairs - coordinates - associated with that bin. Store a list of these lists.
    bin_coordinates = []
    for n in range(nr):
        coord_arrays = np.where(binned_separations==n)
        row_coords = coord_arrays[0]
        col_coords = coord_arrays[1]
        bin_coordinates += [list(zip(row_coords,col_coords))]

    #Divide up the bins if desired. A larger number of jobs is helpful for judging progress as the program is running.
    #N_bin_splits = N_processes//nr + 1
    if N_bin_splits > 1:
        split_size = N_skewers//N_bin_splits
        splits = []
        for split_n in range(N_bin_splits):
            splits += [[split_size*split_n,split_size*(split_n+1)-1]]
        splits[N_bin_splits-1][1] = N_skewers-1
        tasks = []
        for bin_n in range(nr):
            for split in splits:
                tasks += [(bin_n,bin_coordinates[bin_n],skewer_rows[:,split[0]:split[1]])]
    else:
        split_size = N_skewers
        tasks = [(bin_n,bin_coordinates[bin_n],skewer_rows) for bin_n in range(nr)]

    print('divided job up into {} tasks, with ~{} skewers in each one.'.format(len(tasks),split_size))

    #Define the worker function.
    def get_xi(bin_n,bin_coordinates,skewer_rows):

        del_squared_chunk = 0
        N_contributions_chunk = 0
        N_skewers = skewer_rows.shape[1]

        for skewer_n in range(N_skewers):

            skewer = skewer_rows[:,skewer_n]

            for coordinates in bin_coordinates:
                i = coordinates[0]
                j = coordinates[1]

                del_squared_chunk += skewer[i]*skewer[j]
                N_contributions_chunk += 1

        if N_contributions_chunk > 0:
            xi_chunk = del_squared_chunk/N_contributions_chunk
        else:
            xi_chunk = 0

        return [bin_n,N_contributions_chunk,xi_chunk,del_squared_chunk]

    #Define a progress-tracking function.
    def log_result(retval):
        bin_n = retval[0]
        N_contributions_bin = retval[1]
        xi_bin = retval[2]
        del_squared_bin = retval[3]

        N_contributions[bin_n] += N_contributions_bin
        del_squared[bin_n] += del_squared_bin

        results.append(retval)
        N_complete = len(results)
        N_tasks = len(tasks)
        N_tasks_digits = int(np.log10(N_tasks)) + 1

        utils.progress_bar(N_complete,N_tasks,start_time)

    #Define an error-tracking function.
    def log_error(retval):
        print('error',retval)

    #Run the multiprocessing pool
    if __name__ == '__main__':
        pool = Pool(processes = N_processes)
        results = []
        N_contributions = np.zeros(nr)
        del_squared = np.zeros(nr)
        start_time = time.time()

        print('calculating xi')
        for task in tasks:
            pool.apply_async(get_xi,task,callback=log_result,error_callback=log_error)

        #xi = pool.starmap(get_bin_xi,tasks)
        #print('making xi took {}s'.format(time.time()-start_time))

        pool.close()
        pool.join()

    print(' ')

    #Combine results and estimate errors
    print('estimating errors')
    xi = np.zeros(nr)
    mean_xi_squared = np.zeros(nr)
    var_xi = np.zeros(nr)

    for bin_n in range(nr):
        weight_sum = 0
        for result in results:
            if result[0] == bin_n:

                result_N_contributions = result[1]

                if result_N_contributions != 0:
                    weight = result_N_contributions/N_contributions[bin_n]
                    result_xi = result[2]
                else:
                    weight = 0
                    result_xi = 0

                weighted_result_xi = weight*result_xi
                xi[bin_n] += weighted_result_xi
                mean_xi_squared[bin_n] += weight*(result_xi**2)

                weight_sum += weight

        var_xi[bin_n] = mean_xi_squared[bin_n] - (xi[bin_n])**2

    return R_binned, xi, var_xi

"""
