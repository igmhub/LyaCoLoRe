import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sciint

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

        print(np.sum(initial_skewer_rows[i,:]),np.sum(final_skewer_rows[i,:]))

    return final_skewer_rows

#
def get_sigma_kms(T_K):

    sigma_kms = 9.1*((T_K/10000)**(0.5))

    return sigma_kms

#
def get_W(x_kms,T_K):

    sigma_kms = get_sigma_kms(T_K)

    W = ((np.sqrt(2*(np.pi)*(sigma_kms**2)))**(-1.)) * np.exp(-(x_kms**2)/(2*sigma_kms**2))

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
def get_tau_0(z,Om=0.3147):

    return tau_0_rows

#
def add_thermal_skewer_RSDs(initial_tau_rows,initial_density_rows,velocity_skewer_rows_dz,z,r_hMpc,max_N_steps=None):

    N_qso = initial_tau_rows.shape[0]
    N_cells = initial_tau_rows.shape[1]

    if max_N_steps == None:
        max_N_steps = N_cells

    final_tau_rows = np.zeros(initial_tau_rows.shape)

    dkms_dhMpc = get_dkms_dhMpc(z)

    x_kms = r_hMpc * dkms_dhMpc

    velocity_skewer_rows_kms = dz_to_v_kms(velocity_skewer_rows_dz)

    T_K_rows = get_T_K(z,initial_density_rows)

    contribution = np.zeros(x_kms.shape)

    for i in range(N_qso):
        for j in range(N_cells):

            #Integration range method
            x_limit = 5000.0
            n_int = 2500
            x_kms_integration = np.linspace(x_kms[j]-x_limit,x_kms[j]+x_limit,n_int)

            #Linearly inpterpolate the velocity skewers
            v_kms_integration = np.interp(x_kms_integration,x_kms,velocity_skewer_rows_kms[i,:])

            #Nearest neighbout interpolation on velocity skewers
            #v_kms_integration = sciint.interp1d(x_kms,velocity_skewer_rows_kms[i,:],kind='nearest',bounds_error=False,fill_value=0)(x_kms_integration)

            W_argument = x_kms[j]*np.ones(n_int) - x_kms_integration - v_kms_integration
            W = get_W(W_argument,np.interp(x_kms_integration,x_kms,T_K_rows[i,:]))

            """

            #Replicating linear method method
            x_limit = 5000.0
            n_int = 2500
            x_kms_integration = np.linspace(x_kms[j]-x_limit,x_kms[j]+x_limit,n_int)

            W = np.zeros(x_kms_integration.shape)

            #Linearly inpterpolate the velocity skewers
            #v_kms_integration = np.interp(x_kms_integration,x_kms,velocity_skewer_rows_kms[i,:])

            #Nearest neighbout interpolation on velocity skewers
            #v_kms_integration = sciint.interp1d(x_kms,velocity_skewer_rows_kms[i,:],kind='nearest',bounds_error=False,fill_value=(velocity_skewer_rows_kms[i,0],velocity_skewer_rows_kms[i,-1]))(x_kms_integration)

            #"delta function" velocities
            v_kms_integration = np.zeros(x_kms_integration.shape)
            k_values = [k for k in range(N_cells) if x_kms[k] > x_kms[j]-x_limit and x_kms[k] < x_kms[j]+x_limit]
            for k in k_values:
                index = np.searchsorted(x_kms_integration,x_kms[k])
                if abs(x_kms[k] - x_kms_integration[index]) <= abs(x_kms[k] - x_kms_integration[index-1]):
                    v_kms_integration[index] = velocity_skewer_rows_kms[i,k]
                else:
                    v_kms_integration[index-1] = velocity_skewer_rows_kms[i,k]

            xpv_kms_integration = x_kms_integration + v_kms_integration

            if j>0:
                W += np.maximum(np.zeros(x_kms_integration.shape),(xpv_kms_integration - x_kms[j-1]*np.ones(n_int))/(x_kms[j]*np.ones(n_int) - x_kms[j-1]*np.ones(n_int))) * (xpv_kms_integration <= x_kms[j])
            else:
                W += np.ones(x_kms_integration.shape) * (xpv_kms_integration <= x_kms[j])

            if j<N_cells-1:
                W += np.maximum(np.zeros(x_kms_integration.shape),(x_kms[j+1]*np.ones(n_int) - xpv_kms_integration))/(x_kms[j+1]*np.ones(n_int) - x_kms[j]*np.ones(n_int)) * (xpv_kms_integration > x_kms[j])
            else:
                W += np.ones(x_kms_integration.shape) * (xpv_kms_integration >= x_kms[j])

            if j>0:
                if j<N_cells-1:
                    W /= 1.#(0.5*(x_kms[j+1] - x_kms[j-1]))

            """


            #Caluclate tau_R afresh
            #tau_R = get_tau_R(x_kms_integration)

            #Linearly inpterpolate the old tau
            tau_R = np.interp(x_kms_integration,x_kms,initial_tau_rows[i,:])

            #Nearest neighbout interpolation on old tau
            #tau_R = sciint.interp1d(x_kms,initial_tau_rows[i,:],kind='nearest',bounds_error=False,fill_value=0)(x_kms_integration)

            #"delta function" tau
            """
            tau_R = np.zeros(x_kms_integration.shape)
            dx = 2*x_limit/n_int
            k_values = [k for k in range(N_cells) if x_kms[k] > x_kms[j]-x_limit and x_kms[k] < x_kms[j]+x_limit]
            for k in k_values:
                index = np.searchsorted(x_kms_integration,x_kms[k])
                if abs(x_kms[k] - x_kms_integration[index]) <= abs(x_kms[k] - x_kms_integration[index-1]):
                    tau_R[index] = initial_tau_rows[i,k] / dx
                else:
                    tau_R[index-1] = initial_tau_rows[i,k] / dx
            """

            final_tau_rows[i,j] = np.trapz(tau_R*W,x_kms_integration)

    return final_tau_rows
