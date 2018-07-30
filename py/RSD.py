import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sciint
import math

import general

lya = 1215.67

#Function to add linear RSDs from the velocity skewers.
#delete?
def add_linear_skewer_RSDs(initial_skewer,velocity_skewer_dz,z):

    N_qso = initial_skewer.shape[0]
    N_cells = initial_skewer.shape[1]

    final_skewer = np.zeros(initial_skewer.shape)

    for i in range(N_qso):
        for j in range(N_cells):
            #Add the dz from the velocity skewers to get a 'new_z' for each cell
            z_cell = z[j]
            dz_cell = velocity_skewer_dz[i,j]
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

            final_skewer[i,j_upper] += w_upper*initial_skewer[i,j]
            final_skewer[i,j_lower] += w_lower*initial_skewer[i,j]

        #print(np.sum(initial_skewer[i,:]),np.sum(final_skewer[i,:]))

    return final_skewer

#
def get_sigma_kms(T_K):

    sigma_kms = 9.1*((T_K/10000)**(0.5))

    return sigma_kms

#
def get_T_K(z,density):

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

    N_cells = density.shape[1]
    T_K = np.zeros(density.shape)

    for j in range(N_cells):
        T_0_value = np.interp(z[j],z_values,T_0_values)
        gm1_value = np.interp(z[j],z_values,gm1_values)

        T_K[:,j] = T_0_value*((density[:,j])**gm1_value)

    return T_K

def S(x):
    S = x*math.erf(x) + 1./(np.sqrt(np.pi))*np.exp(-(x**2))
    return S

def J(x,a,sigma):
    #Technically should have  -(1/(4*a))*(K(a,sigma)-K(-a,sigma)) too.
    #These cancel out, but if we want to include uneven cell size then they won't do.
    b = 1./(np.sqrt(2)*sigma)
    J = (1./(4*a*b))*(S((x+a)*b)-S((x-a)*b))
    return J

#
def add_skewer_RSDs(initial_tau,initial_density,velocity_skewer_dz,z,r_hMpc,thermal=False):

    N_qso = initial_tau.shape[0]
    N_cells = initial_tau.shape[1]
    final_tau = np.zeros(initial_tau.shape)

    #Convert radial distance to a velocity.
    dkms_dhMpc = general.get_dkms_dhMpc(z)
    x_kms = r_hMpc * dkms_dhMpc

    #Calculate the temperature in every cell if we want to include thermal effects.
    if thermal == True:
        T_K = get_T_K(z,initial_density)

    #count = np.zeros(100)
    #total = 0

    for i in range(N_qso):
        for j in range(N_cells):
            #Add the dz from the velocity skewers to get a 'new_z' for each cell
            z_cell = z[j]
            dz_cell = velocity_skewer_dz[i,j]
            new_z_cell = z_cell + dz_cell

            x_kms_cell = x_kms[j]
            if j > 0 and j < N_cells-1:
                cell_size = (x_kms[j+1] - x_kms[j-1])/2.
            elif j == 0:
                cell_size = (x_kms[j+1] - x_kms[j])
            elif j == N_cells:
                cell_size = (x_kms[j] - x_kms[j-1])

            new_r_hMpc = np.interp(new_z_cell,z,r_hMpc)
            new_x_kms_cell = new_r_hMpc * general.get_dkms_dhMpc(new_z_cell)
            #new_x_kms_cell = x_kms_cell

            j_upper = np.searchsorted(x_kms,new_x_kms_cell)
            j_lower = j_upper - 1

            #If we want to include thermal effects, we include contributions to all cells within a chosen x_kms range.
            if thermal == True:
                tau = initial_tau[i,j]

                #Calculate the dispersion as a result of thermal effects.
                sigma_kms = get_sigma_kms(T_K[i,j])

                #Define the x_kms range over which we will add contributions.
                x_kms_rad = cell_size/2. + 5.*sigma_kms
                x_upper_limit = new_x_kms_cell + x_kms_rad
                x_lower_limit = new_x_kms_cell - x_kms_rad

                #If at least one limit is within the skewer, determine which cells we determine weights for.
                if x_upper_limit > x_kms[0] and x_lower_limit < x_kms[-1]:
                    j_upper_limit = general.NN_sorted(x_kms,x_upper_limit)
                    j_lower_limit = general.NN_sorted(x_kms,x_lower_limit)
                    j_values = np.array(list(range(j_lower_limit,j_upper_limit+1)))
                else:
                    j_values = np.array([])

                #total += 1
                #if len(j_values)<100:
                #    count[len(j_values)] += 1

                cell_weights = np.zeros(len(j_values))

                #For each such cell, find the bottom and top of the cell.
                for j_value in j_values:
                    if j_value > 0 and j_value < N_cells-1:
                        bot = (x_kms[j_value-1] + x_kms[j_value])/2. - new_x_kms_cell
                        top = (x_kms[j_value] + x_kms[j_value+1])/2. - new_x_kms_cell
                    elif j_value == 0:
                        bot = x_kms[j_value] - (x_kms[j_value+1] - x_kms[j_value])/2. - new_x_kms_cell
                        top = (x_kms[j_value] + x_kms[j_value+1])/2. - new_x_kms_cell
                    elif j_value == N_cells-1:
                        bot = (x_kms[j_value-1] + x_kms[j_value])/2. - new_x_kms_cell
                        top = x_kms[j_value] + (x_kms[j_value] - x_kms[j_value-1])/2. - new_x_kms_cell

                    #Use the J function to integrate the weight in this range.
                    # TODO: update this to take into account that consecutive cells may not be precisely the same size
                    #weight = (0.5)*(math.erf((np.sqrt(2))*sigma_kms*top) - math.erf((np.sqrt(2))*sigma_kms*bot))
                    weight = J(top,cell_size/2.,sigma_kms) - J(bot,cell_size/2.,sigma_kms)

                    #cell_weights[j_value-j_values[0]] += weight
                    final_tau[i,j_value] += weight*tau

            #If we do not want to include thermal effects, we only allocate to the cell above and the cell below.
            else:
                tau = initial_tau[i,j]
                #If it has moved off the low-z end of the skewer, lower weight is 0
                #Only include an upper weight if it is within 1 cell's width.
                if j_lower < 0:
                    w_lower = 0.
                    if abs(x_kms[0] - new_x_kms_cell) < abs(x_kms[1] - x_kms[0]):
                        w_upper = 1. - abs(x_kms[0] - new_x_kms_cell)/abs(x_kms[1] - x_kms[0])
                        final_tau[i,j_upper] += w_upper*tau
                    else:
                        w_upper = 0.

                #If it has moved off the high-z end of the skewer, upper weight is 0
                #Only include a lower weight if it is within 1 cell's width.
                elif j_upper >= N_cells:
                    w_upper = 0.
                    if abs(x_kms[-1] - new_x_kms_cell) < abs(x_kms[-1] - x_kms[-2]):
                        w_lower = 1. - abs(x_kms[-1] - new_x_kms_cell)/abs(x_kms[-1] - x_kms[-2])
                        final_tau[i,j_lower] += w_lower*tau
                    else:
                        w_lower = 0.

                #Otherwise, split the contribution between the new neighbours, distance weighted.
                else:
                    x_kms_upper = x_kms[j_upper]
                    x_kms_lower = x_kms[j_lower]

                    w_upper = abs(new_x_kms_cell - x_kms_lower)/(x_kms_upper - x_kms_lower)
                    w_lower = abs(new_x_kms_cell - x_kms_upper)/(x_kms_upper - x_kms_lower)

                    final_tau[i,j_upper] += w_upper*tau
                    final_tau[i,j_lower] += w_lower*tau

            #final_tau[i,j_values] += cell_weights*initial_tau[i,j]

    return final_tau


