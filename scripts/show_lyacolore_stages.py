import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d

from lyacolore import utils

lya = utils.lya_rest

basedir = "/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v9/v9.0.0_all_files/"
N_side = 16
pixel = 0
#i_skewer = 11
i_skewers = list(range(24,100))
fontsize = 16
figsize = (12,15)
dpi = 80

#Set style options everywhere.
#plt.rc('text', usetex=True)
plt.rc('font', size=fontsize)

#Plot variables.
stages_1 = ['picca-gaussian', 'picca-density', 'picca-tau-notnorm', 'picca-flux-notnorm']
label_1 =  [r'$\delta_G$',r'1+$\delta$',r'$\tau$',r'$F$']
stages_2 = ['picca-gaussian-colorecell', None, 'picca-tau-noRSD-notnorm', None]
label_2 =  [r'$\delta_C$',None,r'$\tau_\mathrm{noRSD}$',None]
axis_label = ['Gaussian\nfield','Lognormal\ndensity','Optical\ndepth','Transmitted\nflux fraction']
symmetrical = [True, False, False, False]
add_one = [False,True,False,False]
h_lines = {0: [0], 1: [0], 2: [0], 3: [0,1]}
plot_types = ['skewer']
lambda_min = 3800. #Angstroms
lambda_max = 3900. #Angstroms

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

for i_skewer in i_skewers:

    print(i_skewer)
    #Make the figure.
    fig, axs = plt.subplots(N_stages, N_types, sharex=True, figsize=figsize, dpi=dpi, facecolor='w', edgecolor='k')

    for i in range(N_stages):
        dirname = utils.get_dir_name(basedir,pixel)
        filename = utils.get_file_name(dirname,stages_1[i],N_side,pixel,compressed=True)
        h = fits.open(filename)
        if i == 0:
            mockid_1 = h[3].data['THING_ID'][i_skewer]
        else:
            try:
                i_skewer = np.where(h[3].data['THING_ID']==mockid_1)[0][0]
            except:
                print('Skewer',mockid_1,'not found')
                pass
        lambdas = 10**h[2].data
        skewer = h[0].data[:,i_skewer]
        if add_one[i]:
            skewer += 1.
        print(mockid_1)
        axs[i].plot(lambdas,skewer,label=label_1[i],color=style_dict[stages_1[i]]['c'],linestyle=style_dict[stages_1[i]]['ls'])
        if stages_2[i]:
            filename = utils.get_file_name(dirname,stages_2[i],N_side,pixel,compressed=True)
            h_2 = fits.open(filename)
            mockid_2 = h[3].data['THING_ID'][i_skewer]
            lambdas_2 = 10**h_2[2].data
            skewer_2 = h_2[0].data[:,i_skewer]
            axs[i].plot(lambdas_2,skewer_2,label=label_2[i],color=style_dict[stages_2[i]]['c'],linestyle=style_dict[stages_2[i]]['ls'])
            print(mockid_2)

            if 'RSD' in stages_2[i]:
                filename = utils.get_file_name(dirname,'gaussian-colore',N_side,pixel,compressed=True)
                h = fits.open(filename)
                lambdas_vel = lya*(1+h[4].data['Z'][:])
                i_col = np.where(h[1].data['MOCKID']==mockid_2)[0][0]
                #print(i_col)
                vel = h[3].data[i_col,:]
                axs[i].plot(lambdas_vel,vel*5000,label='vel*5000',color='grey',linestyle=style_dict[stages_2[i]]['ls'])

        print(' ')

        #Scale the axes nicely given the x axis limits.
        j_min = np.searchsorted(lambdas,lambda_min)
        j_max = np.searchsorted(lambdas,lambda_max)
        lambda_range = lambda_max - lambda_min
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
        axs[i].set_ylabel(axis_label[i],rotation=90)
        axs[i].yaxis.set_label_coords(-0.1, 0.5)

        #Add a stage number/section reference:
        #axs[i].text(lambda_min+lambda_range*0.05, y_max-y_range*0.05, 'Stage {}'.format(i), bbox={'facecolor': 'gray', 'alpha': 1.0, 'pad': 4}, verticalalignment='top', horizontalalignment='left')

        #Add a grid and labels.
        #axs[i].grid()
        axs[i].legend(loc=1)

        #Add an arrow to show progression
        if i<N_stages-1:
             axs[i].annotate('', xy=(-0.13, -0.1), xycoords='axes fraction', xytext=(-0.13, 0.1), arrowprops=dict(arrowstyle="->", color='k'))

        #Add horizontal lines
        for h_val in h_lines[i]:
            axs[i].axhline(y=h_val,color='gray',zorder=0,alpha=0.5)

    plt.xlim(lambda_min,lambda_max)
    plt.xlabel(r'$\lambda\ [\mathrm{\AA}]$')
    #fig.align_ylabels()

    m = fits.open(basedir+'/master.fits')
    l = (1+m['COSMO_EXP'].data['Z'])*utils.lya_rest
    r = m['COSMO_EXP'].data['R']
    r_to_l = interp1d(l,r)
    l_to_r = interp1d(r,l)
    topax = axs[0].twiny()

    """
    print(l[:5])
    print(r[:5])
    print(' ')
    print(l[5:10])
    print(r[5:10])
    print(' ')
    """

    topax_values = [3940,3980,4020,4060]
    topax_locations = r_to_l(topax_values)

    #topax_locations = axs[0].get_xticks()
    #topax_values = l_to_r(topax_locations)

    topax_strs = [str(round(xval,0)) for xval in topax_values]
    topax.set_xlim(axs[0].get_xlim())
    topax.set_xticks(topax_locations)
    topax.set_xticklabels(topax_strs)
    topax.set_xlabel(r'$r\ [\rm{Mpc}/h]$')

    fig.tight_layout()
    fig.subplots_adjust(hspace=0)
    plt.savefig('skewers.pdf')
    plt.show()
