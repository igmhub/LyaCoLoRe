from . import absorber

"""
This is a dictionary containing names of the added metals, their corresponding restframe wavelength, and optical depth.
The values are taken as the ones used in desisim, and Si lines were chosen as in Bautista et al. 2017 https://arxiv.org/pdf/1702.00176.pdf
Notice that the metals that are being considered are only those with wavelengths  close from Lya.

"""
def get_lyb_absorber():
    return absorber.AbsorberData('Lyb', 1025.72, 0.1901)

def get_metal_dict(selection=None,metals_list=None):

    if (selection is None) and (metals_list is None):
        selection = 'standard'

    if metals_list is None:
        if selection=='standard':
            metals_list = ['CIV(1551)', 'CIV(1548)']
            #print(metals_list)
        elif selection=='full':
            metals_list = list(metals_data.keys())

    metals_dict = {}
    for m in metals_list:
        metals_dict[m] = absorber.AbsorberData(m,
                metals_data[m]['rest_wave'],metals_data[m]['strength'],
                metals_data[m]['HDU_name'])

    return metals_dict

metals_data = {
    ## The below metals form the "standard" selection.
    'SiII(1260)' : {'rest_wave': 1260.4221000, 'strength': 3.54200E-04, 'HDU_name': 'SI1260'},
    'SiIII(1207)': {'rest_wave': 1206.5000000, 'strength': 1.89190E-03, 'HDU_name': 'SI1207'},
    'SiII(1193)' : {'rest_wave': 1193.2897000, 'strength': 9.07760E-04, 'HDU_name': 'SI1193'},
    'SiII(1190)' : {'rest_wave': 1190.4158000, 'strength': 1.28478E-03, 'HDU_name': 'SI1190'},

    'MgI(2853)'  : {'rest_wave': 2852.9600000, 'strength': 1.00000E-04, 'HDU_name': 'MG2853'},
    'MgII(2804)' : {'rest_wave': 2803.5324000, 'strength': 5.00000E-04, 'HDU_name': 'MG2804'},
    'MgII(2796)' : {'rest_wave': 2796.3511000, 'strength': 9.00000E-04, 'HDU_name': 'MG2796'},
    'FeII(2600)' : {'rest_wave': 2600.1724835, 'strength': 1.00000E-04, 'HDU_name': 'FE2600'},
    'FeII(2587)' : {'rest_wave': 2586.6495659, 'strength': 1.00000E-04, 'HDU_name': 'FE2587'},
    'MnII(2577)' : {'rest_wave': 2576.8770000, 'strength': 1.00000E-04, 'HDU_name': 'MN2577'},
    'FeII(2383)' : {'rest_wave': 2382.7641781, 'strength': 1.00000E-04, 'HDU_name': 'FE2383'},
    'FeII(2374)' : {'rest_wave': 2374.4603294, 'strength': 1.00000E-04, 'HDU_name': 'FE2374'},
    'FeII(2344)' : {'rest_wave': 2344.2129601, 'strength': 1.00000E-04, 'HDU_name': 'FE2344'},
    'AlIII(1863)': {'rest_wave': 1862.7911300, 'strength': 1.00000E-04, 'HDU_name': 'AL1863'},
    'AlIII(1855)': {'rest_wave': 1854.7182900, 'strength': 1.00000E-04, 'HDU_name': 'AL1855'},
    'AlII(1671)' : {'rest_wave': 1670.7886000, 'strength': 1.00000E-04, 'HDU_name': 'AL1671'},
    'FeII(1608)' : {'rest_wave': 1608.4511000, 'strength': 1.00000E-04, 'HDU_name': 'FE1608'},
    'CIV(1551)'  : {'rest_wave': 1550.7784500, 'strength': 5.43500E-04, 'HDU_name': 'C1551' },
    'CIV(1548)'  : {'rest_wave': 1548.2049000, 'strength': 1.48700E-03, 'HDU_name': 'C1548' },
    'SiII(1527)' : {'rest_wave': 1526.7069800, 'strength': 1.00000E-04, 'HDU_name': 'SI1527'},
    'SiIV(1403)' : {'rest_wave': 1402.7729100, 'strength': 5.00000E-04, 'HDU_name': 'SI1403'},
    'SiIV(1394)' : {'rest_wave': 1393.7601800, 'strength': 9.00000E-04, 'HDU_name': 'SI1394'},
    'CII(1335)'  : {'rest_wave': 1334.5323000, 'strength': 1.00000E-04, 'HDU_name': 'C1335' },
    'SiII(1304)' : {'rest_wave': 1304.3702000, 'strength': 1.00000E-04, 'HDU_name': 'SI1304'},
    'OI(1302)'   : {'rest_wave': 1302.1685000, 'strength': 1.00000E-04, 'HDU_name': 'O1302' },
    'NV(1243)'   : {'rest_wave': 1242.8040000, 'strength': 5.00000E-04, 'HDU_name': 'N1243' },
    'NV(1239)'   : {'rest_wave': 1238.8210000, 'strength': 5.00000E-04, 'HDU_name': 'N1239' },
    'NI(1200)'   : {'rest_wave': 1200.0000000, 'strength': 1.00000E-03, 'HDU_name': 'N1200' },
    'OI(1039)'   : {'rest_wave': 1039.2300000, 'strength': 1.00000E-03, 'HDU_name': 'O1039' },
    'OVI(1038)'  : {'rest_wave': 1037.6130000, 'strength': 3.82000E-01, 'HDU_name': 'O1038' },
    'OVI(1032)'  : {'rest_wave': 1031.9120000, 'strength': 5.35800E-03, 'HDU_name': 'O1032' },
    'CIII(977)'  : {'rest_wave':  977.0200000, 'strength': 5.00000E-03, 'HDU_name': 'C977'  },
    'OI(989)'    : {'rest_wave':  988.7000000, 'strength': 1.00000E-03, 'HDU_name': 'O989'  },
    'SiII(990)'  : {'rest_wave':  989.8731000, 'strength': 1.00000E-03, 'HDU_name': 'SI990' },
    'LY3'        : {'rest_wave':  972.5370000, 'strength': 6.97000E-02, 'HDU_name': 'LY3'   },
    'LY4'        : {'rest_wave':  949.7431000, 'strength': 3.35000E-02, 'HDU_name': 'LY4'   },
    'LY5'        : {'rest_wave':  937.8035000, 'strength': 1.87000E-02, 'HDU_name': 'LY5'   },
    }
