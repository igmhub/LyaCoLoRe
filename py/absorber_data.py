from . import absorber

"""
This is a dictionary containing names of the added metals, their corresponding restframe wavelength, and optical depth.
The values are taken as the ones used in desisim, and Si lines were chosen as in Bautista et al. 2017 https://arxiv.org/pdf/1702.00176.pdf
Notice that the metals that are being considered are only those with wavelengths  close from Lya.

"""
def get_lyb_absorber():
    return absorber.AbsorberData('Lyb', 1025.72, 0.1901)

def get_metal_dict():
    metal_dict = {
      'SiII(1260)'  : absorber.AbsorberData('SiII(1260)',  1260.42, 3.542e-4),
      'SiIII(1207)' : absorber.AbsorberData('SiIII(1207)', 1206.50, 1.8919e-3),
      'SiII(1193)'  : absorber.AbsorberData('SiII(1193)',  1193.29, 9.0776e-4),
      'SiII(1190)'  : absorber.AbsorberData('SiII(1190)',  1190.42, 6.4239e-4),
      #Until the memory needs are known, the following metals will not be added.
      #'NV(1243)'    : absorber.AbsorberData('NV(1243)',  1242.804, 5.e-4),
      #'NV(1239)'    : absorber.AbsorberData('NV(1239)',  1238.821, 5.e-4),
      #'NI(1200)'    : absorber.AbsorberData('NI(1200)',  1200.,    1.e-3),
      #'OI(1039)'    : absorber.AbsorberData('OI(1039)',  1039.230, 1.e-3),
      #'OVI(1038)'   : absorber.AbsorberData('OVI(1038)', 1037.613, 3.382-3),
      #'OVI(1032)'   : absorber.AbsorberData('OVI(1032)', 1031.912, 5.358e-3),
    }
    return metal_dict
