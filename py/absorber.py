import numpy as np

from . import utils

# Class storing information related to a particular absorber
class AbsorberData:

    #Initialisation function.
    def __init__(self,name='LYA',rest_wave=utils.lya_rest,flux_transform_m=1.0):
        self.name = name
        self.rest_wave = rest_wave
        self.flux_transform_m = flux_transform_m

        # we will store here the optical depth and record when RSDs have been applied
        self.tau = None
        self.RSDs_applied = False

        return

    def tau_computed(self):
        if self.tau is None: 
            return False
        return True

    def transmission(self):
        if not self.tau_computed():
            print('you can not get transmission without first computing tau')
            raise ValueError('Can not access transmission arrays')
        return np.exp(-self.tau)
