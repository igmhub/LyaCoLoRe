# Class storing information related to a particular absorber
class AbsorberData:
    #Initialisation function.
    def __init__(self,name='LYA',rest_wave=1215.67,flux_transform_m=1.0):
        self.name = name
        self.rest_wave = rest_wave
        self.flux_transform_m = flux_transform_m
        return
