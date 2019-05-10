import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

from lyacolore import utils

basedir = "/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v6/v6.0.0/" 
N_side = 16
pixel = 0
i_skewer = 0

#Plot variables.
stages_1 = ['picca-gaussian', 'picca-density', 'picca-tau-notnorm', 'picca-flux-notnorm']
label_1 = [r'$\delta_G$',r'1+$\delta$',r'$\tau$',r'$F$']
stages_2 = ['picca-gaussian-colorecell', None, 'picca-tau-noRSD-notnorm', 'picca-flux-noRSD-notnorm']
label_2 = [r'$\delta_C$',None,r'$\tau_\mathrm{noRSD}$',r'$F_\mathrm{noRSD}$']

symmetrical = [True, False, False, False]
add_one = [False,True,False,False]
plot_types = ['skewer']
lambda_min = 3750. #Angstroms
lambda_max = 3850. #Angstroms

style_dict = {'picca-gaussian-colorecell': {'c': 'C0', 'ls': '--'},
              'picca-gaussian':            {'c': 'C0', 'ls': '-'},
              'picca-density':             {'c': 'C1', 'ls': '-'},
              'picca-tau-noRSD-notnorm':   {'c': 'C2', 'ls': ':'},
              'picca-tau-notnorm':         {'c': 'C2', 'ls': '-'},
              'picca-flux-noRSD-notnorm':  {'c': 'C3', 'ls': ':'},
              'picca-flux-notnorm':        {'c': 'C3', 'ls': '-'},
              }

#Deduced variables.
N_stages = len(stages_1)
N_types = len(plot_types)

#Make the subplots, and reduce the horizontal space between axes to 0.
fig, axs = plt.subplots(N_stages, N_types, sharex=True, figsize=(8, 10), dpi= 80, facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace=0)

for i in range(N_stages):
    dirname = utils.get_dir_name(basedir,pixel)
    filename = utils.get_file_name(dirname,stages_1[i],N_side,pixel)
    h = fits.open(filename)
    lambdas = 10**h[2].data
    skewer = h[0].data[:,i_skewer]
    if add_one[i]:
        skewer += 1.
    print(h[3].data['THING_ID'][i_skewer])
    axs[i].plot(lambdas,skewer,label=label_1[i],color=style_dict[stages_1[i]]['c'],linestyle=style_dict[stages_1[i]]['ls'])
    if stages_2[i]:
        filename = utils.get_file_name(dirname,stages_2[i],N_side,pixel)
        h_2 = fits.open(filename)
        lambdas_2 = 10**h_2[2].data
        skewer_2 = h_2[0].data[:,i_skewer]
        axs[i].plot(lambdas_2,skewer_2,label=label_2[i],color=style_dict[stages_2[i]]['c'],linestyle=style_dict[stages_2[i]]['ls'])
        print(h_2[3].data['THING_ID'][i_skewer])
    print(' ')

    #Scale the axes nicely given the x axis limits.
    j_min = np.searchsorted(lambdas,lambda_min)
    j_max = np.searchsorted(lambdas,lambda_max)
    y_min = np.min(skewer[j_min:j_max])
    y_max = np.max(skewer[j_min:j_max])
    y_range = y_max - y_min
    y_low_lim = y_min - y_range*0.1
    y_upp_lim = y_max + y_range*0.1
    if symmetrical[i]:
        y_low_lim = -np.max((abs(y_low_lim),abs(y_upp_lim)))
        y_upp_lim = np.max((abs(y_low_lim),abs(y_upp_lim)))
    axs[i].set_ylim(y_low_lim, y_upp_lim)
    axs[i].tick_params(which='both')

    #Add a grid and labels.
    axs[i].grid()
    axs[i].legend(loc=1)

plt.xlim(lambda_min,lambda_max)
plt.xlabel(r'$\lambda\ /\ \mathrm{\AA}$')
plt.savefig('skewers.pdf')
plt.show()
