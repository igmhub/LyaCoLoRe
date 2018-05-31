import numpy as np

def beta_QSO_kaiser(z,b,Om_z0=0.3147):

    Om = Om_z0 * ((1+z)**3) / (Om_z0 * ((1+z)**3) + 1 - Om_z0)
    Ol = 1 - Om
    # TODO: replace this with the analytic expression using D from file
    f = (Om**0.6) + (Ol/70.)*(1 + Om/2.) #using review eqn 4.14
    beta = f/b

    return beta
