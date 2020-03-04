from . import absorber

"""
This is a dictionary containing names of the added metals, their corresponding restframe wavelength, and optical depth.
The values are taken as the ones used in desisim, and Si lines were chosen as in Bautista et al. 2017 https://arxiv.org/pdf/1702.00176.pdf
Notice that the metals that are being considered are only those with wavelengths  close from Lya.

"""
def get_lyb_absorber():
    return absorber.AbsorberData('Lyb', 1025.72, 0.1901, 'LYB')

def get_metal_dict():
    metal_dict = {

    # OLD DATA
    #Used in version 7.x. Boosted by 50% in v8.x.x following analysis.
    #'SiII(1190)'  : absorber.AbsorberData('SiII(1190)',  1190.42, 6.4239e-4, 'SI1190'),

    'SiII(1260)'  : absorber.AbsorberData('SiII(1260)', 1260.4221,    3.54200e-04, 'SI1260'),
    'SiIII(1207)' : absorber.AbsorberData('SiIII(1207)',1206.5,       1.89190e-03, 'SI1207'),
    'SiII(1193)'  : absorber.AbsorberData('SiII(1193)', 1193.2897,    9.07760e-04, 'SI1193'),
    'SiII(1190)'  : absorber.AbsorberData('SiII(1190)', 1190.4158,    1.28478e-03, 'SI1190'),

    'MgI(2853)'   : absorber.AbsorberData('MgI(2853)',  2852.96,      1.00000e-04, 'MG2853'),
    'MgII(2804)'  : absorber.AbsorberData('MgII(2804)', 2803.5324,    5.00000e-04, 'MG2804'),
    'MgII(2796)'  : absorber.AbsorberData('MgII(2796)', 2796.3511,    9.00000e-04, 'MG2796'),
    'FeII(2600)'  : absorber.AbsorberData('FeII(2600)', 2600.1724835, 1.00000e-04, 'FE2600'),
    'FeII(2587)'  : absorber.AbsorberData('FeII(2587)', 2586.6495659, 1.00000e-04, 'FE2587'),
    'MnII(2577)'  : absorber.AbsorberData('MnII(2577)', 2576.877,     1.00000e-04, 'MN2577'),
    'FeII(2383)'  : absorber.AbsorberData('FeII(2383)', 2382.7641781, 1.00000e-04, 'FE2383'),
    'FeII(2374)'  : absorber.AbsorberData('FeII(2374)', 2374.4603294, 1.00000e-04, 'FE2374'),
    'FeII(2344)'  : absorber.AbsorberData('FeII(2344)', 2344.2129601, 1.00000e-04, 'FE2344'),
    'AlIII(1863)' : absorber.AbsorberData('AlIII(1863)',1862.79113,   1.00000e-04, 'AL1863'),
    'AlIII(1855)' : absorber.AbsorberData('AlIII(1855)',1854.71829,   1.00000e-04, 'AL1855'),
    'AlII(1671)'  : absorber.AbsorberData('AlII(1671)', 1670.7886,    1.00000e-04, 'AL1671'),
    'FeII(1608)'  : absorber.AbsorberData('FeII(1608)', 1608.4511,    1.00000e-04, 'FE1608'),
    'CIV(1551)'   : absorber.AbsorberData('CIV(1551)',  1550.77845,   5.43500e-04, 'C1551'),
    'CIV(1548)'   : absorber.AbsorberData('CIV(1548)',  1548.2049,    1.48700e-03, 'C1548'),
    'SiII(1527)'  : absorber.AbsorberData('SiII(1527)', 1526.70698,   1.00000e-04, 'SI1527'),
    'SiIV(1403)'  : absorber.AbsorberData('SiIV(1403)', 1402.77291,   5.00000e-04, 'SI1403'),
    'SiIV(1394)'  : absorber.AbsorberData('SiIV(1394)', 1393.76018,   9.00000e-04, 'SI1394'),
    'CII(1335)'   : absorber.AbsorberData('CII(1335)',  1334.5323,    1.00000e-04, 'C1335'),
    'SiII(1304)'  : absorber.AbsorberData('SiII(1304)', 1304.3702,    1.00000e-04, 'SI1304'),
    'OI(1302)'    : absorber.AbsorberData('OI(1302)',   1302.1685,    1.00000e-04, 'O1302'),
    'NV(1243)'    : absorber.AbsorberData('NV(1243)',   1242.804,     5.00000e-04, 'N1243'),
    'NV(1239)'    : absorber.AbsorberData('NV(1239)',   1238.821,     5.00000e-04, 'N1239'),
    'NI(1200)'    : absorber.AbsorberData('NI(1200)',   1200.0,       1.00000e-03, 'N1200'),
    'OI(1039)'    : absorber.AbsorberData('OI(1039)',   1039.23,      1.00000e-03, 'O1039'),
    'OVI(1038)'   : absorber.AbsorberData('OVI(1038)',  1037.613,     3.82000e-01, 'O1038'),
    'OVI(1032)'   : absorber.AbsorberData('OVI(1032)',  1031.912,     5.35800e-03, 'O1032'),
    'CIII(977)'   : absorber.AbsorberData('CIII(977)',  977.02,       5.00000e-03, 'C977'),
    'OI(989)'     : absorber.AbsorberData('OI(989)',    988.7,        1.00000e-03, 'O989'),
    'SiII(990)'   : absorber.AbsorberData('SiII(990)',  989.8731,     1.00000e-03, 'SI990'),
    'LY3'         : absorber.AbsorberData('LY3',        972.537,      6.97000e-02, 'LY3'),
    'LY4'         : absorber.AbsorberData('LY4',        949.7431,     3.35000e-02, 'LY4'),
    'LY5'         : absorber.AbsorberData('LY5',        937.8035,     1.87000e-02, 'LY5'),
    }

    return metal_dict
