import absorber

def get_metal_dict():
    metal_dict = { 
      'SiII(1260)' : absorber.AbsorberData('SiII(1260)',1260.42,8.e-4),
      'SiIII(1207)' : absorber.AbsorberData('SiIII(1207)',1206.50,5.e-3),
      'SiII(1193)' : absorber.AbsorberData('SiII(1193)',1193.29,5.e-4),
      'SiII(1190)' : absorber.AbsorberData('SiII(1190)',1190.42,5.e-4)
    }
    return metal_dict
