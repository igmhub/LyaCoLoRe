import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import lya_mock_functions as mock
import sys
import itertools
import process_functions as functions

delta_types = ['gaussian','density','flux']
base_dir ='/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_revamp/process_output_4096_32/'
h = fits.open(base_dir + 'nside_8_master.fits')
pixel_list = list(sorted(set(h[1].data['PIXNUM'])))
pixel_list = list(range(250))
h.close()

N_side = 8

def get_mean(DELTA_columns,IVAR_columns):
    WEIGHTED_DELTA_columns = DELTA_columns*IVAR_columns
    N_skewers = DELTA_columns.shape[1]
    SUM_column = np.sum(WEIGHTED_DELTA_columns,axis=1)
    N_column = np.sum(IVAR_columns,axis=1)
    MEAN_column = SUM_column/N_column
    return N_skewers,MEAN_column

"""
tasks = []
for delta_type in delta_types:
    for pixel in pixels:
        tasks += [(delta_type,pixel)]
"""

for dt,delta_type in enumerate(delta_types):

    print('looking at {}'.format(delta_type))

    fig = plt.figure()
    N_skewers_total = 0
    max_skewer_length = 0
    picca_mean_total = np.zeros(max_skewer_length)
    picca_lambdas_total = np.zeros(0)

    for p,pixel in enumerate(pixel_list):

        print('looking at pixel {}'.format(pixel),end='\r')

        pixel_100 = int(pixel/100)
        picca = fits.open('{}/{}/{}/picca-{}-{}-{}.fits'.format(base_dir,str(pixel_100),str(pixel),delta_type,N_side,pixel))
        picca_deltas = picca[0].data
        picca_ivars = picca[1].data
        picca_lambdas = 10**picca[2].data
        picca.close()

        N_skewers,picca_mean = get_mean(picca_deltas,picca_ivars)
        plt.plot(picca_lambdas,picca_mean,color='0.9')

        skewer_length_deficiency = max_skewer_length - picca_mean.shape[0]
        max_skewer_length = max(max_skewer_length,picca_mean.shape[0])

        if skewer_length_deficiency > 0:
            picca_mean = np.concatenate((picca_mean,np.zeros(skewer_length_deficiency)))
        elif skewer_length_deficiency < 0:
            picca_mean_total = np.concatenate((picca_mean_total,np.zeros(int(abs(skewer_length_deficiency)))))
            picca_lambdas_total = picca_lambdas

        picca_mean_total = (picca_mean_total*N_skewers_total + picca_mean*N_skewers)/(N_skewers_total + N_skewers)

        N_skewers_total += N_skewers

    plt.plot(picca_lambdas_total,picca_mean_total)
    plt.plot(picca_lambdas_total,np.zeros(picca_lambdas_total.shape[0]),color='0.0')
    plt.plot(picca_lambdas_total,np.average(picca_mean_total)*np.ones(picca_lambdas_total.shape[0]),color='0.0')

    upper_limit = max(abs(picca_mean_total))*1.5
    lower_limit = min(min(picca_mean_total),-upper_limit/3.)*1.5
    plt.ylim(lower_limit,upper_limit)

    plt.title(delta_type)
    plt.savefig('{}_{}_{}.pdf'.format(delta_type,pixel_list[0],pixel_list[-1]))

    print(' ')


"""
print('looking at physical colore')

fig = plt.figure()
N_skewers_total = 0
max_skewer_length = 0
colore_mean_total = np.zeros(max_skewer_length)
colore_lambdas_total = np.zeros(0)

for p,pixel in enumerate(pixel_list):

    print('looking at pixel {}'.format(pixel),end='\r')

    pixel_100 = int(pixel/100)
    colore = fits.open('{}/{}/{}/physical-colore-{}-{}.fits'.format(base_dir,str(pixel_100),str(pixel),N_side,pixel))
    colore_deltas = picca[0].data
    #NEED TO MAKE ALL THE S~TUFF FOR COLORE IVARS
    colore_ivars = functions.make_IVAR_rows(1215.67,Z_QSO,LOGLAM_MAP,N_qso,N_cells)
    colore_lambdas = 10**picca[2].data
    colore.close()

    N_skewers,colore_mean = get_mean(colore_deltas,colore_ivars)
    plt.plot(colore_lambdas,colore_mean,color='0.9')

    skewer_length_deficiency = max_skewer_length - picca_mean.shape[0]
    max_skewer_length = max(max_skewer_length,picca_mean.shape[0])

    if skewer_length_deficiency > 0:
        picca_mean = np.concatenate((picca_mean,np.zeros(skewer_length_deficiency)))
    elif skewer_length_deficiency < 0:
        colore_mean_total = np.concatenate((picca_mean_total,np.zeros(int(abs(skewer_length_deficiency)))))
        colore_lambdas_total = colore_lambdas

    colore_mean_total = (colore_mean_total*N_skewers_total + picca_mean*N_skewers)/(N_skewers_total + N_skewers)

    N_skewers_total += N_skewers


plt.plot(picca_lambdas_total,picca_mean_total)
plt.plot(picca_lambdas_total,np.zeros(picca_lambdas_total.shape[0]),color='0.0')
plt.plot(picca_lambdas_total,np.average(picca_mean_total)*np.ones(picca_lambdas_total.shape[0]),color='0.0')

limit = max(abs(picca_mean_total))*1.5
plt.ylim(-limit,limit)

plt.title(delta_type)
"""

plt.show()
