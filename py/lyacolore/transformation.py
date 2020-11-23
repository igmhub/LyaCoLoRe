import numpy as np
import json

def randomise_parameter_values(data,seed=0):

    print('INFO: adding random noise to non-fixed initial parameters.')

    # Set the seed.
    np.random.seed(seed=seed)

    # For each tuning quantity:
    for k in data.keys():

        # For each parameter defining that quantity:
        for p in data[k]['parameters'].keys():
            p_dict = data[k]['parameters'][p]

            # If desired, also add limits.
            if ('lower_limit' in p_dict.keys()) and ('upper_limit' in p_dict.keys()):
                limit = (p_dict['lower_limit'], p_dict['upper_limit'])
            elif ('lower_limit' in p_dict.keys()):
                limit = (p_dict['lower_limit'], float("infinity"))
            elif ('upper_limit' in p_dict.keys()):
                limit = (-float("infinity"), p_dict['upper_limit'])
            else:
                limit = (-float("infinity"), float("infinity"))

            if (not p_dict['fix']):
                lo = (limit[0]-p_dict['value'])/p_dict['error']
                hi = (limit[1]-p_dict['value'])/p_dict['error']
                loc = p_dict['value']
                scale = p_dict['error']
                p_dict['value'] = truncnorm.rvs(lo,hi,loc=loc,scale=scale,size=1)[0]

    return data

# TODO: At the moment, the Transformation object is really just a set of functions. It needs to contain the information about how to actually carry out the transformations. I think it should be merged with convert, and maybe independent too.

#Object to store the transformation used to go from CoLoRe's Gaussian skewers to
#flux skewers. Stores functions of redshift.
class Transformation:
    def __init__(self):
        return

    @classmethod
    def make_transformation_from_tuning_file(cls,filepath,random=False,seed=0):

        # Initiate the transformation object.
        transformation = cls()

        # Open the tuning file.
        with open(filepath, 'r') as json_file:
            data = json.load(json_file)

        # Randomise values if desired.
        if random:
            randomise_parameter_values(data,seed=seed)

        # Make a dictionary of quantities, for which we have a
        # TransformationQuantity for each key.
        quantities = {}
        for key in data.keys():
            parameters = {k:data[key]['parameters'][k]['value'] for k in data[key]['parameters'].keys()}
            tuning_quantity = TransformationQuantity(key,data[key]['functiontype'],parameters)
            quantities[key] = tuning_quantity

        transformation.quantities = quantities

        return transformation

    def update_parameters_from_minuit(self,minuit_values):

        # For each quantity in this transformation:
        for q in self.quantities.keys():

            # Find the minuit values which correspond to parameters for this
            # quantity.
            relevant_values = [k for k in minuit_values.keys() if k[:len(q)]==q]

            # Make a {parameter name: parameter value} dict for the parameter
            # relevant to this quantity.
            parameters = {p.split('-')[-1]: minuit_values[p] for p in relevant_values}

            # Update the quantity to use the new parameters.
            self.quantities[q].update(parameters)

        return


class TransformationQuantity:
    def __init__(self,name,function_type,parameters):
        self.name = name
        self.function_type = function_type
        self.parameters = parameters
        self.parameter_history = [parameters]
        self._function = function_types[self.function_type](**self.parameters)
        return

    def __call__(self,x):
        return self._function(x)

    def update(self,new_parameters):
        self.parameters = new_parameters
        self.parameter_history += [new_parameters]
        self._function = function_types[self.function_type](**self.parameters)
        return


class Constant:
    def __init__(self,A0=1.):
        self.A0 = A0
        return
    def __call__(self,x):
        if not np.isscalar(x):
            return self.A0*np.ones_like(x)
        else:
            return self.A0

class QuadraticLog:
    def __init__(self,A0=1.,A1=0.,A2=0.,z0=3.):
        self.A0 = A0
        self.A1 = A1
        self.A2 = A2
        self.z0 = z0
        return
    def __call__(self,z):
        x = (1+z)/(1+self.z0)
        val = np.log(self.A0) + self.A1*np.log(x) + self.A2*(np.log(x))**2
        return np.exp(val)

class BasicPowerFunction:
    def __init__(self,n=0.,k1_kms=1.):
        self.n = n
        self.k1_kms = k1_kms
        return
    def __call__(self,k_kms):
        # power used to make mocks in from McDonald et al. (2006)
        val = ((1.0+pow(0.01/self.k1_kms,self.n)) / (1.0+pow(k_kms/self.k1_kms,self.n)))
        return val

function_types = {'basic_power':     BasicPowerFunction,
                  'quadratic_log':   QuadraticLog,
                  'constant':        Constant,
}
