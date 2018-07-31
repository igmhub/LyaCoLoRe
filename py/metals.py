import absorber



"""
This is a dictionary containing names of the added metals, their corresponding restframe wavelength, and optical depth.
The values are taken as the ones used in desisim, and Si lines were chosen as in Bautista et al. 2017 https://arxiv.org/pdf/1702.00176.pdf
Notice that the metals that are being considered are only those with wavelengths  close from Lya.

"""
def get_metal_dict():
    metal_dict = { 
      'SiII(1260)' : absorber.AbsorberData('SiII(1260)',1260.42,8.e-4),
      'SiIII(1207)' : absorber.AbsorberData('SiIII(1207)',1206.50,5.e-3),
      'SiII(1193)' : absorber.AbsorberData('SiII(1193)',1193.29,5.e-4),
      'SiII(1190)' : absorber.AbsorberData('SiII(1190)',1190.42,5.e-4)
    }
    return metal_dict
