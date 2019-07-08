import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.sparse import dok_matrix,csc_matrix,csr_matrix
import math
import time

from . import utils

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

#
def S(x):
    S = x*math.erf(x) + 1./(np.sqrt(np.pi))*np.exp(-(x**2))
    return S

#
def J(x,a,sigma):
    #Technically should have  -(1/(4*a))*(K(a,sigma)-K(-a,sigma)) too.
    #These cancel out, but if we want to include uneven cell size then they won't do.
    b = 1./(np.sqrt(2)*sigma)
    J = (1./(4*a*b))*(S((x+a)*b)-S((x-a)*b))
    return J

#Integral of the erf function in a finite range.
def int_erf(upp,low):
    int = S(upp) - S(low)
    return int

#Weight assigned to an interval (a,b) given a starting velocity xstar, Gaussian
#velocity dispersion sigma and cell width 2d.
def W_interval(a,b,xstar,sigma,d):

    upp1 = (b - xstar + d)/(np.sqrt(2) * sigma)
    low1 = (a - xstar + d)/(np.sqrt(2) * sigma)
    upp2 = (b - xstar - d)/(np.sqrt(2) * sigma)
    low2 = (a - xstar - d)/(np.sqrt(2) * sigma)

    W = (sigma/(2*np.sqrt(2)*d)) * (int_erf(upp1,low1) - int_erf(upp2,low2))

    return W

#
def add_skewer_RSDs(initial_tau,initial_density,velocity_skewer_dz,z,r_hMpc,z_qso,thermal=False,weights=None,d=0.0,z_r0=2.5):

    N_qso = initial_tau.shape[0]
    N_cells = initial_tau.shape[1]

    #Convert radial distance to a velocity.
    dkms_dhMpc = utils.get_dkms_dhMpc(z)
    x_kms = r_hMpc * dkms_dhMpc

    #Calculate the temperature in every cell if we want to include thermal effects.
    if thermal == True:
        T_K = get_T_K(z,initial_density)

    #Get the RSD weights if we don't already have them.
    if not weights:
        weights = get_weights(initial_density,velocity_skewer_dz,z,r_hMpc,z_qso,thermal=thermal,d=d,z_r0=z_r0)

    #Compute the final values of tau.
    start = time.time()

    final_tau = np.zeros(initial_tau.shape)
    for k in range(N_qso):
        skewer_weights = weights[k]
        final_tau[k,:] = skewer_weights.dot(initial_tau[k,:].T)

    return final_tau

#
def get_weights(initial_density,velocity_skewer_dz,z,r_hMpc,z_qso,thermal=False,d=0.0,z_r0=2.5):

    #Get the pivot value of r.
    r0 = np.interp(z_r0,z,r_hMpc)

    #Set up variables.
    N_qso = velocity_skewer_dz.shape[0]
    N_cells = velocity_skewer_dz.shape[1]
    weights = {}

    #Convert radial distance to a velocity.
    dkms_dhMpc = utils.get_dkms_dhMpc(z)
    x_kms = r_hMpc * dkms_dhMpc

    #Calculate the temperature in every cell if we want to include thermal effects.
    if thermal == True:
        T_K = get_T_K(z,initial_density)

    #Calculate the edges of the cells in terms of r, z and x.
    r_edges = utils.get_edges(r_hMpc)
    z_edges = interp1d(r_hMpc,z,fill_value='extrapolate')(r_edges)
    x_edges = interp1d(r_hMpc,x_kms,fill_value='extrapolate')(r_edges)

    #Define the upper and lower edges in terms of x, and get the cell sizes.
    x_uedges = x_edges[1:]
    x_ledges = x_edges[:-1]
    old_cell_sizes = x_uedges - x_ledges

    start = time.time()

    for i in range(N_qso):

        #Set up the variables for storing the sparse matrix data.
        indices = []
        data = []
        indptr = [0]

        #extract the velocity skewer, and shift the edges of the cells.
        dz = velocity_skewer_dz[i,:]

        #Calculate dz of the edges by getting midpoints of dz of cells.
        #This takes into account the compression/expansion of the cells.
        #dz_edges_mid = utils.get_edges(dz)
        #dz_ledges = dz_edges_mid[:-1]
        #dz_uedges = dz_edges_mid[1:]

        #Calculated dz of the edges from the dz of the cells.
        #This ignores the compression/expansion of the cells.
        dz_ledges = dz
        dz_uedges = dz

        #Shift all of the z edges and cells.
        z_ledges_shifted = z_edges[:-1] + dz_ledges
        z_uedges_shifted = z_edges[1:] + dz_uedges
        z_shifted = z + dz

        #Calculate the shifted r edges and cells by interpolating.
        r_of_z = interp1d(z_edges,r_edges,fill_value='extrapolate')
        r_ledges_shifted = r_of_z(z_ledges_shifted)
        r_uedges_shifted = r_of_z(z_uedges_shifted)
        r_shifted = r_of_z(z_shifted)

        #Shift the cell again to simulate an extra velocity gradient.
        r_ledges_shifted -= (r_hMpc - r0) * d
        r_uedges_shifted -= (r_hMpc - r0) * d
        r_shifted -= (r_hMpc - r0) * d

        #Calculate the shifted x edges and cells by interpolating.
        x_of_r = interp1d(r_edges,x_edges,fill_value='extrapolate')
        x_ledges_shifted = x_of_r(r_ledges_shifted)
        x_uedges_shifted = x_of_r(r_uedges_shifted)
        x_shifted = x_of_r(r_shifted)

        #Get the new cell sizes.
        new_cell_sizes = x_uedges_shifted - x_ledges_shifted

        #Go through each cell up to the cell the QSO is in.
        j_limit = np.searchsorted(z_edges[:-1],z_qso[i])
        for j in range(j_limit):

            #Define new variables for commonly used vector elements.
            new_x_kms_cell_ledge = x_ledges_shifted[j]
            new_x_kms_cell_uedge = x_uedges_shifted[j]
            new_x_kms_cell = x_shifted[j]
            old_cell_size = old_cell_sizes[j]
            new_cell_size = new_cell_sizes[j]

            #If we want to include thermal effects, we include contributions to all cells within a chosen x_kms range.
            if thermal == True:

                #Calculate the dispersion as a result of thermal effects.
                sigma_kms = get_sigma_kms(T_K[i,j])

                #Define the x_kms range over which we will add contributions.
                x_kms_rad = old_cell_size/2. + 5.*sigma_kms
                x_upper_limit = new_x_kms_cell + x_kms_rad
                x_lower_limit = new_x_kms_cell - x_kms_rad

                #If at least one limit is within the skewer, determine which cells we determine weights for.
                j_values = np.where((x_upper_limit > x_kms) * (x_lower_limit < x_kms))[0]

                #For each such cell, find the bottom and top of the cell.
                #Shift so that the centre of the new cell is at 0.
                for k,j_value in enumerate(j_values):

                    """
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
                    weight_old = J(top,cell_size/2.,sigma_kms) - J(bot,cell_size/2.,sigma_kms)
                    """

                    weight = W_interval(x_ledges[j_value],x_uedges[j_value],new_x_kms_cell,sigma_kms,cell_size/2.)
                    #weight = J(x_uedges[j_value]-new_x_kms_cell,cell_size/2.,sigma_kms) - J(x_ledges[j_value]-new_x_kms_cell,cell_size/2.,sigma_kms)

                    indices += [j_value]
                    data += [weight]

                indptr += [(indptr[-1]+len(j_values))]

            #If we do not want to include thermal effects, we only allocate to
            #the cell above and the cell below according to the overlap.
            else:

                #If the cell ends up having some overlap with the skewer, find which cells it contributes to.
                j_values = list(np.where((new_x_kms_cell_uedge>x_ledges)*(new_x_kms_cell_ledge<x_uedges))[0])

                #For the cells it contributes to, find how its weight is spread
                #cell_weights = np.maximum(0.,np.minimum(x_uedges[j_values],new_x_kms_cell_uedge) - np.maximum(x_ledges[j_values],new_x_kms_cell_ledge)) / (new_x_kms_cell_uedge - new_x_kms_cell_ledge)
                #w = list(cell_weights)

                w = []
                for j_value in j_values:
                    overlap = max(0., min(x_uedges[j_value], new_x_kms_cell_uedge) - max(x_ledges[j_value], new_x_kms_cell_ledge))
                    weight = overlap/new_cell_size
                    w += [weight]

                #Add the data to the inputs for the sparse matrix
                indices += j_values
                indptr += [(indptr[-1] + len(j_values))]
                data += w

        #Make the sparse matrix, and add it to the dictionary.
        indptr += [indptr[-1]]*(N_cells + 1 - len(indptr))
        csc_weights = csc_matrix((data, indices, indptr), shape=(N_cells,N_cells))
        weights[i] = csc_weights

    #print('time to calc RSDs weights is {:2.3f}'.format(time.time()-start))

    return weights


################################################################################

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
