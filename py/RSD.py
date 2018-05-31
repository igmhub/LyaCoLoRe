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

#
def get_sigma_kms(T_K):

    sigma_kms = 9.1*((T_K/10000)**(0.5))

    return sigma_kms

#
def get_W(x_kms,T_K):

    sigma_kms = get_sigma_kms(T_K)

    W = (np.sqrt(2*(np.pi)*(sigma_kms**2))) * np.exp(-(x_kms**2)/(2*sigma_kms**2))

    return W

#
def get_T_K(z,density_rows):

    #Data from McDonald et al. 2001
    #Probably more up to date versions available
    T_14_McD = [20700.,20300.,20100.]
    T_0_McD =  [17400.,18400.,17400.]
    gm1_McD =  [0.52  ,0.29  ,0.43  ]
    z_McD =    [2.4   ,3.0   ,3.9   ]

    #Data from Hiss et al. 2017
    T_0_Hiss = [13530.,9580.,14272.,17454.,23130.,20608.,19172.,17348.]
    gm1_Hiss = [0.39  ,0.59 ,0.49  ,0.37  ,0.31  ,0.13  ,0.48  ,0.42  ]
    z_Hiss =   [2.0   ,2.2  ,2.4   ,2.6   ,2.8   ,3.0   ,3.2   ,3.4   ]

    #Data from Ricotti et al. 2000
    T_0_Ric =  [4700. ,17700.,25200.,19800.]
    gm1_Ric =  [0.85  ,0.32  ,0.22  ,0.38  ]
    z_Ric =    [0.06  ,1.9   ,2.75  ,3.56  ]

    #Mixed data
    T_0_values = [17700.,13530.,9580.,14272.,17454.,23130.,20608.,19172.,17348.,19800.,20100.]
    gm1_values = [0.32  ,0.39  ,0.59 ,0.49  ,0.37  ,0.31  ,0.13  ,0.48  ,0.42  ,0.38  ,0.43  ]
    z_values =   [1.9   ,2.0   ,2.2  ,2.4   ,2.6   ,2.8   ,3.0   ,3.2   ,3.4   ,3.56  ,3.9   ]

    N_cells = density_rows.shape[1]
    T_K = np.zeros(density_rows.shape)

    for j in range(N_cells):
        T_0_value = np.interp(z[j],z_values,T_0_values)
        gm1_value = np.interp(z[j],z_values,gm1_values)

        T_K[:,j] = T_0_value*((density_rows[:,j])**gm1_value)

    return T_K

#This should probably go in a general function file
def get_dkms_dhMpc(z,Om=0.3147):

    E_z = np.sqrt(Om*(1+z)**3 + (1-Om))
    dkms_dhMpc = 100. * E_z / (1+z)

    return dkms_dhMpc

#Converts a peculiar redshift to a peculiar velocity in kms-1.
def dz_to_v_kms(dz):

    c = 299792.458 #kms-1
    alpha = (1 + dz)**2
    v_kms = c*((alpha-1)/(alpha+1))

    return v_kms

#
def get_tau_real(vel_real_kms,tau_redshift,vel_redshift_kms,v_parallel_kms,T_K):

    W_integrand = np.zeros(tau_redshift.shape)

    for i in range(tau_redshift.shape[0]):
        W_integrand[i] = get_W(vel_real_kms-vel_redshift_kms[i]-v_parallel_kms[i],T_K[i])

    tau_real = np.trapz(tau_redshift*W_integrand,vel_redshift_kms)

    return tau_real

#
def add_thermal_RSDs(initial_tau_rows,initial_density_rows,velocity_skewer_rows_dz,z,r_hMpc,max_N_steps=None):

    N_qso = initial_tau_rows.shape[0]
    N_cells = initial_tau_rows.shape[1]

    if max_N_steps == None:
        max_N_steps = N_cells

    final_tau_rows = np.zeros(initial_tau_rows.shape)

    dkms_dhMpc = get_dkms_dhMpc(z)

    x_kms = r_hMpc * dkms_dhMpc

    velocity_skewer_rows_kms = dz_to_v_kms(velocity_skewer_rows_dz)

    T_K_rows = get_T_K(z,initial_density_rows)

    for i in range(N_qso):
        for j in range(N_cells):

            #Given the max number of steps to integrate over, work out which j values we will integrate over.
            #Could also do this by specifying a range of r or x or z?
            j_values = list(range(max(0,j-max_N_steps),min(j+max_N_steps+1,N_cells)))
            x_kms_integration = x_kms[j_values]

            W_argument = x_kms[j]*np.ones(len(j_values)) - x_kms_integration - velocity_skewer_rows_kms[i,j_values]

            W = get_W(W_argument,T_K_rows[i,j_values])

            #I think we can just use the tau we have already here?
            #tau_R = get_tau_R(x_kms_integration)
            tau_R = initial_tau_rows[i,j_values]

            final_tau_rows[i,j] = np.trapz(tau_R*W,x_kms_integration)

    return final_tau_rows
