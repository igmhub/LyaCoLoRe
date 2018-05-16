import numpy as np

#Function to add linear RSDs from the velocity skewers.
def add_linear_skewer_RSDs(initial_skewer_rows,velocity_skewer_rows_dz,z):

    N_qso = initial_skewer_rows.shape[0]
    N_cells = initial_skewer_rows.shape[1]

    final_skewer_rows = np.zeros(initial_skewer_rows.shape)

    for i in range(N_qso):
        for j in range(N_cells):
            #Add the dz from the velocity skewers to get a 'new_z' for each cell
            z_cell = z[j]
            dz_cell = velocity_skewer_rows_dz[i,j]
            new_z_cell = z_cell + dz_cell

            #Work out where in the skewer the cell 'moves' to.
            #i.e. what are the new neighbouring cells.
            j_upper = np.searchsorted(z,new_z_cell)
            j_lower = j_upper - 1

            #If it has moved off the z=0 end of the skewer, push it back.
            if j_lower < 0:
                w_upper = 1.0
                w_lower = 0.0
                j_lower += 1

            #If it has moved off the max z end of the skewer, push it back.
            elif j_upper >= N_cells:
                w_lower = 1.0
                w_upper = 0.0
                j_upper -= 1

            #Otherwise, split the contribution between the new neighbours, distance weighted.
            else:
                z_upper = z[j_upper]
                z_lower = z[j_lower]

                w_upper = abs(new_z_cell - z_lower)/(z_upper - z_lower)
                w_lower = abs(new_z_cell - z_upper)/(z_upper - z_lower)

            final_skewer_rows[i,j_upper] += w_upper*initial_skewer_rows[i,j]
            final_skewer_rows[i,j_lower] += w_lower*initial_skewer_rows[i,j]

    return final_skewer_rows


def get_sigma_kms(T_K):

    sigma_kms = 9.1*((T_K/10000)**(0.5))

    return sigma_kms


def get_W(x_kms,T_K):

    sigma_kms = get_sigma_kms(T_K)

    W = (np.sqrt(2*(np.pi)*(sigma_kms**2))) * np.exp(-(x_kms**2)/(2*sigma_kms**2))

    return W


def get_T_K(z,rho):

    #Data from McDonald et al. 2001
    #Probably more up to date versions available
    T_14_McD = [20700.,20300.,20100.]
    gm1_McD = [0.52,0.29,0.43]
    z_McD = [2.4,3.0,3.9]

    T_14_values = np.interp(z,z_McD,T_14_McD)
    gm1_values = np.interp(z,z_McD,gm1_McD)

    T_K = T_14_value*((rho/1.4)**gm1_values)

    return T_K

#This should probably go in a general function file
def get_dkms_dhMpc(z,Om=0.3147):

    E_z = np.sqrt(Om*(1+z)**3 + (1-Om))
    dkms_dhMpc = 100. * E_z / (1+z)

    return dkms_dhMpc


def get_tau_real(vel_real_kms,tau_redshift,vel_redshift_kms,v_parallel_kms,T_K):

    W_integrand = np.zeros(tau_redshift.shape)

    for i in range(tau_redshift.shape[0]):
        W_integrand[i] = get_W(vel_real_kms-vel_redshift_kms[i]-v_parallel_kms[i],T_K[i])

    tau_real = np.trapz(tau_redshift*W_integrand,vel_redshift_kms)

    return tau_real


def add_thermal_RSDs(initial_tau_rows,initial_density_rows,velocity_skewer_rows_dz,z,x_hMpc,neighbours=False):

    N_cells = initial_tau_rows.shape[0]
    N_qso = initial_tau_rows.shape[1]

    final_tau_rows = np.zeros(initial_tau_rows)

    dkms_dhMpc = get_dkms_dhMpc(z)

    for i in range(N_qso):

        T_K = get_T_K(z,initial_density_rows[i,:])

        for j in range(N_cells):

            #Do we alraedy have a cosmological velocity? There's a V column in colore output but not sure what it is
            x_kms = x_hMpc[j]*dkms_dhMpc[j]

            if neighbours == True:
                js = list(range(max(0,j-1),min(j+1,N_cells)))
                final_tau_rows[i,j] = get_tau_real(x_kms,initial_tau_rows[i,js],x_hMpc[js]*dkms_dhMpc[js],velocity_skewer_rows_kms[js],T_K)
            else:
                final_tau_rows[i,j] = get_tau_real(x_kms,initial_tau_rows[i,:],x_hMpc[:]*dkms_dhMpc[:],velocity_skewer_rows_kms[:],T_K)

    return
