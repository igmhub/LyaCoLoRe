import numpy as np

#Function to define the lognormal transform.
def lognormal_transform(y,sG,D):
    return np.exp(D*y - 0.5*(D**2)*(sG**2))

#Function to convert gaussian field skewers (in rows) to lognormal delta skewers (in rows).
def gaussian_to_lognormal_delta(GAUSSIAN_DELTA_rows,SIGMA_G,D):

    LN_DENSITY_rows = np.zeros(GAUSSIAN_DELTA_rows.shape)
    SIGMA_G = SIGMA_G*np.ones(GAUSSIAN_DELTA_rows.shape[1])
    D = D*np.ones(GAUSSIAN_DELTA_rows.shape[1])

    for j in range(GAUSSIAN_DELTA_rows.shape[1]):
        LN_DENSITY_rows[:,j] = np.exp(D[j]*GAUSSIAN_DELTA_rows[:,j] - ((D[j])**2)*(SIGMA_G[j]**2)/2.)#lognormal_transform(GAUSSIAN_DELTA_rows[:,j],SIGMA_G[j],D[j])

    LN_DENSITY_DELTA_rows = LN_DENSITY_rows - 1.

    #Not sure why this was there?
    #LN_DENSITY_DELTA_rows = LN_DENSITY_DELTA_rows.astype('float32')

    return LN_DENSITY_DELTA_rows

#Function to convert from density to tau using alpha*density^beta
def density_to_tau(density,alpha,beta):

    tau = alpha*(density**beta)

    return tau

#Function to convert from density to flux using e^-(alpha*density^beta)
def density_to_flux(density,alpha,beta):

    tau = density_to_tau(density,alpha,beta)
    F = np.exp(-tau)

    return F

#Function to convert lognormal delta skewers (in rows) to gaussian field skewers (in rows).
def lognormal_delta_to_gaussian(LN_DENSITY_DELTA_rows,SIGMA_G,D):

    LN_DENSITY_rows = 1.0 + LN_DENSITY_DELTA_rows

    GAUSSIAN_DELTA_rows = np.zeros(LN_DENSITY_DELTA_rows.shape)

    if np.array(D).shape[0] == 1:
        D = np.ones(LN_DENSITY_rows.shape[1])*D

    if np.array(SIGMA_G).shape[0] == 1:
        SIGMA_G = np.ones(LN_DENSITY_rows.shape[1])*SIGMA_G

    for j in range(GAUSSIAN_DELTA_rows.shape[1]):
        GAUSSIAN_DELTA_rows[:,j] = (np.log(LN_DENSITY_rows[:,j]))/D[j] + (D[j])*(SIGMA_G[j]**2)/2

    GAUSSIAN_DELTA_rows = GAUSSIAN_DELTA_rows.astype('float32')

    return GAUSSIAN_DELTA_rows

#Function to convert tau skewers (in rows) to density skewers (in rows).
def tau_to_density(TAU_rows,alpha,beta):

    DENSITY_rows = np.zeros(TAU_rows.shape)
    N_cells = TAU_rows.shape[1]

    for j in range(N_cells):
        DENSITY_rows[:,j] = (TAU_rows[:,j]/alpha[j])**(1/beta)

    return DENSITY_rows
