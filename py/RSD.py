import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sciint
import math

import general

lya = 1215.67


#Function to add linear RSDs from the velocity skewers.
#delete?
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

        #print(np.sum(initial_skewer_rows[i,:]),np.sum(final_skewer_rows[i,:]))

    return final_skewer_rows

#
def get_sigma_kms(T_K):

    sigma_kms = 9.1*((T_K/10000)**(0.5))

    return sigma_kms

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
def add_skewer_RSDs(initial_tau_rows,initial_density_rows,velocity_skewer_rows_dz,z,r_hMpc,thermal=False):

    N_qso = initial_tau_rows.shape[0]
    N_cells = initial_tau_rows.shape[1]
    final_tau_rows = np.zeros(initial_tau_rows.shape)

    #Convert radial distance to a velocity.
    dkms_dhMpc = get_dkms_dhMpc(z)
    x_kms = r_hMpc * dkms_dhMpc

    #Calculate the temperature in every cell if we want to include thermal effects.
    if thermal == True:
        T_K_rows = get_T_K(z,initial_density_rows)

    #count = np.zeros(100)
    #total = 0

    for i in range(N_qso):
        for j in range(N_cells):
            #Add the dz from the velocity skewers to get a 'new_z' for each cell
            z_cell = z[j]
            dz_cell = velocity_skewer_rows_dz[i,j]
            new_z_cell = z_cell + dz_cell

            x_kms_cell = x_kms[j]
            if j > 0 and j < N_cells-1:
                cell_size = (x_kms[j+1] - x_kms[j-1])/2.
            elif j == 0:
                cell_size = (x_kms[j+1] - x_kms[j])
            elif j == N_cells:
                cell_size = (x_kms[j] - x_kms[j-1])

            new_r_hMpc = np.interp(new_z_cell,z,r_hMpc)
            new_x_kms_cell = new_r_hMpc * get_dkms_dhMpc(new_z_cell)
            #new_x_kms_cell = x_kms_cell

            j_upper = np.searchsorted(x_kms,new_x_kms_cell)
            j_lower = j_upper - 1

            #If we want to include thermal effects, we include contributions to all cells within a chosen x_kms range.
            if thermal == True:
                tau = initial_tau_rows[i,j]

                #Calculate the dispersion as a result of thermal effects.
                sigma_kms = get_sigma_kms(T_K_rows[i,j])

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
                    final_tau_rows[i,j_value] += weight*tau

            #If we do not want to include thermal effects, we only allocate to the cell above and the cell below.
            else:
                tau = initial_tau_rows[i,j]
                #If it has moved off the low-z end of the skewer, lower weight is 0
                #Only include an upper weight if it is within 1 cell's width.
                if j_lower < 0:
                    w_lower = 0.
                    if abs(x_kms[0] - new_x_kms_cell) < abs(x_kms[1] - x_kms[0]):
                        w_upper = 1. - abs(x_kms[0] - new_x_kms_cell)/abs(x_kms[1] - x_kms[0])
                        final_tau_rows[i,j_upper] += w_upper*tau
                    else:
                        w_upper = 0.

                #If it has moved off the high-z end of the skewer, upper weight is 0
                #Only include a lower weight if it is within 1 cell's width.
                elif j_upper >= N_cells:
                    w_upper = 0.
                    if abs(x_kms[-1] - new_x_kms_cell) < abs(x_kms[-1] - x_kms[-2]):
                        w_lower = 1. - abs(x_kms[-1] - new_x_kms_cell)/abs(x_kms[-1] - x_kms[-2])
                        final_tau_rows[i,j_lower] += w_lower*tau
                    else:
                        w_lower = 0.

                #Otherwise, split the contribution between the new neighbours, distance weighted.
                else:
                    x_kms_upper = x_kms[j_upper]
                    x_kms_lower = x_kms[j_lower]

                    w_upper = abs(new_x_kms_cell - x_kms_lower)/(x_kms_upper - x_kms_lower)
                    w_lower = abs(new_x_kms_cell - x_kms_upper)/(x_kms_upper - x_kms_lower)

                    final_tau_rows[i,j_upper] += w_upper*tau
                    final_tau_rows[i,j_lower] += w_lower*tau

            #final_tau_rows[i,j_values] += cell_weights*initial_tau_rows[i,j]

    return final_tau_rows


"""
# TODO: DELETE THE BELOW WHEN TIDYING
BELOW IS OLD CODE.
INCLUDES CODE TO USE AN INTEGRAL METHOD FOR DETERMINING THERMAL RSDS.
SLOWER THAN THE METHOD INCLUDED ABOVE

def I(x,a,sigma,N_steps=10**4):

    x_int = np.linspace(-a,a,N_steps+1)
    # TODO: should it be 1/2N+1 rather than 1/2N?
    integrand = (1/(2*a*np.sqrt(2*np.pi)*sigma)) * np.exp(-((x - x_int)**2)/(2*sigma**2))
    I = np.trapz(integrand,)

    return I

def K(x,sigma):
    K = x*math.erf(x/(np.sqrt(2)*sigma)) + ((np.sqrt(2)*sigma)/(np.sqrt(np.pi)))*np.exp(-(x**2)/(2*(sigma**2)))
    return K



#
def get_W(x_kms,T_K):

    sigma_kms = get_sigma_kms(T_K)

    W = ((np.sqrt(2*(np.pi)*(sigma_kms**2)))**(-1.)) * np.exp(-(x_kms**2)/(2*sigma_kms**2))

    return W


#
def get_tau_0(z,Om=0.3147):

    return tau_0_rows


#
def get_tau_real(vel_real_kms,tau_redshift,vel_redshift_kms,v_parallel_kms,T_K):

    W_integrand = np.zeros(tau_redshift.shape)

    for i in range(tau_redshift.shape[0]):
        W_integrand[i] = get_W(vel_real_kms-vel_redshift_kms[i]-v_parallel_kms[i],T_K[i])

    tau_real = np.trapz(tau_redshift*W_integrand,vel_redshift_kms)

    return tau_real


#Don't think this is correct
#Converts a peculiar redshift to a peculiar velocity in kms-1.
def dz_to_dv_kms(dz,z,v_kms):

    c = 299792.458 #kms-1
    alpha = (1 + z + dz)
    dv_kms = ((alpha**2)*(c - v_kms) - (c + v_kms))/(1 + alpha**2)

    return dv_kms



velocity_skewer_rows_kms = np.zeros(velocity_skewer_rows_dz.shape)
for i in range(N_qso):
    velocity_skewer_rows_kms[i,:] = dz_to_dv_kms(velocity_skewer_rows_dz[i,:],z,x_kms)

linear_cell_weights_kms = np.zeros(len(j_values))

j_upper = np.searchsorted(x_kms,new_x_kms_cell)
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
    x_kms_upper = x_kms[j_upper]
    x_kms_lower = x_kms[j_lower]

    w_upper = abs(new_x_kms_cell - x_kms_lower)/(x_kms_upper - x_kms_lower)
    w_lower = abs(new_x_kms_cell - x_kms_upper)/(x_kms_upper - x_kms_lower)

linear_cell_weights_kms[j_upper - j_values[0]] = w_upper
linear_cell_weights_kms[j_lower - j_values[0]] = w_lower

cell_weights = linear_cell_weights_kms

"""


"""

linear_cell_weights_z = np.zeros(len(j_values))

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

linear_cell_weights_z[j_upper - j_values[0]] = w_upper
linear_cell_weights_z[j_lower - j_values[0]] = w_lower

cell_weights = linear_cell_weights_z

#print(j)

#print(z_cell,dz_cell,new_z_cell)
#print(x_kms_cell,vel_kms_cell,new_x_kms_cell)

#print(z[j_values])
#print(x_kms[j_values])

#print(linear_cell_weights_z)
#print(linear_cell_weights_kms)
#print(cell_weights)

#print(' ')

"""

"""

#Integration Method

#Integration range method
x_limit = x_kms_rad
n_int = 2500
x_kms_integration = np.linspace(x_kms[j]-x_limit,x_kms[j]+x_limit,n_int+1)

x_kms_integration = x_kms_integration[x_kms_integration>0]
x_kms_integration = x_kms_integration[x_kms_integration<=x_kms[-1]]
n_int = x_kms_integration.shape[0]

#Linearly inpterpolate the velocity and temperature skewers
v_kms_integration = np.interp(x_kms_integration,x_kms,velocity_skewer_rows_kms[i,:],left=0.0,right=0.0)
T_K_integration = np.interp(x_kms_integration,x_kms,T_K_rows[i,:])

#Nearest neighbout interpolation on velocity skewers
#v_kms_integration = sciint.interp1d(x_kms,velocity_skewer_rows_kms[i,:],kind='nearest',bounds_error=False,fill_value=0)(x_kms_integration)
#T_K_integration = sciint.interp1d(x_kms,T_K_rows[i,:],kind='nearest',bounds_error=False,fill_value=0)(x_kms_integration)

W_argument = x_kms[j]*np.ones(n_int) - x_kms_integration - v_kms_integration
W = get_W(W_argument,T_K_integration)

#
#if i==6 and j in list(range(marker-2,marker+3)):
#    print(j,x_kms[0],x_kms_integration[0],'//',x_kms[j],'//',x_kms[-1],x_kms_integration[-1])

#    W_argument_nv = x_kms[j]*np.ones(n_int) - x_kms_integration
#    W_nv = get_W(W_argument_nv,np.interp(x_kms_integration,x_kms,T_K_rows[i,:]))
#    A = np.array(list(zip(x_kms_integration,W,W_nv,v_kms_integration)))
#    np.savetxt('/Users/jfarr/Projects/test_data/W_data/W_{}_{}.txt'.format(i,j),A)

#    A = np.array(list(zip(x_kms,velocity_skewer_rows_kms[i,:])))
#    np.savetxt('/Users/jfarr/Projects/test_data/W_data/vel_{}_{}.txt'.format(i,j),A)

#Replicating linear method method
#x_limit = 5000.0
#n_int = 2500
#x_kms_integration = np.linspace(x_kms[j]-x_limit,x_kms[j]+x_limit,n_int)

#W = np.zeros(x_kms_integration.shape)

#Linearly inpterpolate the velocity skewers
#v_kms_integration = np.interp(x_kms_integration,x_kms,velocity_skewer_rows_kms[i,:])

#Nearest neighbout interpolation on velocity skewers
#v_kms_integration = sciint.interp1d(x_kms,velocity_skewer_rows_kms[i,:],kind='nearest',bounds_error=False,fill_value=(velocity_skewer_rows_kms[i,0],velocity_skewer_rows_kms[i,-1]))(x_kms_integration)

#"delta function" velocities
#v_kms_integration = np.zeros(x_kms_integration.shape)
#k_values = [k for k in range(N_cells) if x_kms[k] > x_kms[j]-x_limit and x_kms[k] < x_kms[j]+x_limit]
#for k in k_values:
#    index = np.searchsorted(x_kms_integration,x_kms[k])
#    if abs(x_kms[k] - x_kms_integration[index]) <= abs(x_kms[k] - x_kms_integration[index-1]):
#        v_kms_integration[index] = velocity_skewer_rows_kms[i,k]
#    else:
#        v_kms_integration[index-1] = velocity_skewer_rows_kms[i,k]

#xpv_kms_integration = x_kms_integration + v_kms_integration

#if j>0:
#    W += np.maximum(np.zeros(x_kms_integration.shape),(xpv_kms_integration - x_kms[j-1]*np.ones(n_int))/(x_kms[j]*np.ones(n_int) - x_kms[j-1]*np.ones(n_int))) * (xpv_kms_integration <= x_kms[j])
#else:
#    W += np.ones(x_kms_integration.shape) * (xpv_kms_integration <= x_kms[j])

#if j<N_cells-1:
#    W += np.maximum(np.zeros(x_kms_integration.shape),(x_kms[j+1]*np.ones(n_int) - xpv_kms_integration))/(x_kms[j+1]*np.ones(n_int) - x_kms[j]*np.ones(n_int)) * (xpv_kms_integration > x_kms[j])
#else:
#    W += np.ones(x_kms_integration.shape) * (xpv_kms_integration >= x_kms[j])

#if j>0:
#    if j<N_cells-1:
#        W /= 1.#(0.5*(x_kms[j+1] - x_kms[j-1]))


#Caluclate tau_R afresh
#tau_R = get_tau_R(x_kms_integration)

#Linearly inpterpolate the old tau
tau_R = np.interp(x_kms_integration,x_kms,initial_tau_rows[i,:],left=0.0,right=0.0)

#Nearest neighbout interpolation on old tau
#tau_R = sciint.interp1d(x_kms,initial_tau_rows[i,:],kind='nearest',bounds_error=False,fill_value=0)(x_kms_integration)

#"delta function" tau

#tau_R = np.zeros(x_kms_integration.shape)
#dx = 2*x_limit/n_int
#k_values = [k for k in range(N_cells) if x_kms[k] > x_kms[j]-x_limit and x_kms[k] < x_kms[j]+x_limit]
#for k in k_values:
#    index = np.searchsorted(x_kms_integration,x_kms[k])
#    if abs(x_kms[k] - x_kms_integration[index]) <= abs(x_kms[k] - x_kms_integration[index-1]):
#        tau_R[index] = initial_tau_rows[i,k] / dx
#    else:
#        tau_R[index-1] = initial_tau_rows[i,k] / dx

final_tau_rows[i,j] = np.trapz(tau_R*W,x_kms_integration)
"""

        #print(np.trapz(final_tau_rows[i,:],z)/np.trapz(initial_tau_rows[i,:],z))
    #print('{:3.1%} {:3.1%} {:3.1%} {:3.1%} {:3.1%} {:3.1%} {:3.1%} {:3.1%} {:3.1%} {:3.1%} {:3.1%}'.format(count[0]/total,count[1]/total,count[2]/total,count[3]/total,count[4]/total,count[5]/total,count[6]/total,count[7]/total,count[8]/total,count[9]/total,count[10]/total))
