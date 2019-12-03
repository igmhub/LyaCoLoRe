import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d

from lyacolore import utils

lya = utils.lya_rest

basedir = "/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v9/v9.0.0_all_files/"
N_side = 16
pixel = 0
fontsize = 16
figsize = (12,15)
dpi = 80
y_buffer = 0.1

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

colours = ['#F5793A','#A95AA1','#85C0F9','#0F2080']

style_dict = {'picca-gaussian-colorecell': {'c': colours[0], 'ls': '--'},
              'picca-gaussian':            {'c': colours[0], 'ls': '-'},
              'picca-density':             {'c': colours[1], 'ls': '-'},
              'picca-tau-noRSD-notnorm':   {'c': colours[2], 'ls': ':'},
              'picca-tau-notnorm':         {'c': colours[2], 'ls': '-'},
              'picca-flux-noRSD-notnorm':  {'c': colours[3], 'ls': ':'},
              'picca-flux-notnorm':        {'c': colours[3], 'ls': '-'},
              }

#Deduced variables.
N_stages = len(stages_1)
N_types = len(plot_types)
dirname = utils.get_dir_name(basedir,pixel)

files_1 = []
files_2 = []
for i in range(N_stages):
    filename = utils.get_file_name(dirname,stages_1[i],N_side,pixel,compressed=True)
    h = fits.open(filename)
    files_1 += [h]
    if stages_2[i] is not None:
        filename = utils.get_file_name(dirname,stages_2[i],N_side,pixel,compressed=True)
        h = fits.open(filename)
        files_2 += [h]
    else:
        files_2 += [None]

mockids = files_1[0][3].data['THING_ID']
#mockids = [64347,110595,110963,111104]
mockids = [64347]

def plot_skewer(ax,h,mockid,label,c,ls,add_one=False):
    lambdas = 10**h[2].data
    try:
        i_skewer = np.where(h[3].data['THING_ID']==mockid)[0][0]
        skewer = h[0].data[:,i_skewer]
        if add_one:
            skewer += 1
        ax.plot(lambdas,skewer,label=label,color=c,linestyle=ls)
        return lambdas,skewer
    except:
        print('Skewer',mockid,'not found')
        return lambdas,None


for mockid in mockids:

    print(mockid)

    #Make the figure.
    fig, axs = plt.subplots(N_stages, N_types, sharex=True, figsize=figsize, dpi=dpi, facecolor='w', edgecolor='k')

    for i in range(N_stages):
        lambdas,skewer = plot_skewer(axs[i],files_1[i],mockid,label_1[i],style_dict[stages_1[i]]['c'],style_dict[stages_1[i]]['ls'],add_one=add_one[i])
        if stages_2[i] is not None:
            lambdas_2,skewer_2 = plot_skewer(axs[i],files_2[i],mockid,label_2[i],style_dict[stages_2[i]]['c'],style_dict[stages_2[i]]['ls'],add_one=add_one[i])

        #Scale the axes nicely given the x axis limits.
        j_min = np.searchsorted(lambdas,lambda_min)
        j_max = np.searchsorted(lambdas,lambda_max)
        lambda_range = lambda_max - lambda_min
        if skewer is not None:
            y_min = np.min(skewer[j_min:j_max])
            y_max = np.max(skewer[j_min:j_max])
            y_range = y_max - y_min
            y_low_lim = y_min - y_range*y_buffer
            y_upp_lim = y_max + y_range*y_buffer
            if symmetrical[i]:
                y_low_lim = -np.max((abs(y_low_lim),abs(y_upp_lim)))
                y_upp_lim = np.max((abs(y_low_lim),abs(y_upp_lim)))
            axs[i].set_ylim(y_low_lim, y_upp_lim)
        axs[i].tick_params(which='both')
        axs[i].set_ylabel(axis_label[i],rotation=90)
        axs[i].yaxis.set_label_coords(-0.1, 0.5)

        """
        if stages_2[i] is not None:
            if 'RSD' in stages_2[i]:
                filename = utils.get_file_name(dirname,'gaussian-colore',N_side,pixel,compressed=True)
                h_col = fits.open(filename)
                lambdas_vel = lya*(1+h_col[4].data['Z'][:])
                try:
                    i_col = np.where(h_col[1].data['MOCKID']==mockid)[0][0]
                    vel = h_col[3].data[i_col,:]
                    ax2 = axs[i].twinx()
                    ax2.plot(lambdas_vel,vel,label=r'$\rm{d}z_{\rm{RSD}}$',color='grey',linestyle=style_dict[stages_2[i]]['ls'])
                    j_min_vel = np.searchsorted(lambdas_vel,lambda_min)
                    j_max_vel = np.searchsorted(lambdas_vel,lambda_max)
                    y_min_vel = np.min(vel[j_min_vel:j_max_vel])
                    y_max_vel = np.max(vel[j_min_vel:j_max_vel])
                    y_range_vel = y_max_vel - y_min_vel
                    y_low_lim_vel = y_min_vel - y_range_vel*y_buffer
                    y_upp_lim_vel = y_max_vel + y_range_vel*y_buffer
                    y_low_lim_vel = -np.max((abs(y_low_lim_vel),abs(y_upp_lim_vel)))
                    y_upp_lim_vel = np.max((abs(y_low_lim_vel),abs(y_upp_lim_vel)))
                    y_upp_lim_vel = y_low_lim_vel*y_upp_lim/y_low_lim
                    ax2.set_ylim(y_low_lim_vel,y_upp_lim_vel)
                except IndexError:
                    print('Skewer',mockid,'not found')
        """

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

    #Add a top axis with the distance in Mpc/h
    m = fits.open(basedir+'/master.fits')
    l = (1+m['COSMO_EXP'].data['Z'])*utils.lya_rest
    r = m['COSMO_EXP'].data['R']
    r_to_l = interp1d(l,r)
    l_to_r = interp1d(r,l)
    topax = axs[0].twiny()

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
    plt.savefig('skewers_{}.pdf'.format(mockid))
    plt.show()


"""
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
        try:
            i_skewer_2 = np.where(h[3].data['THING_ID']==mockid_1)[0][0]
        except:
            print('Skewer',mockid_1,'not found')
            pass
        lambdas_2 = 10**h_2[2].data
        skewer_2 = h_2[0].data[:,i_skewer_2]
        axs[i].plot(lambdas_2,skewer_2,label=label_2[i],color=style_dict[stages_2[i]]['c'],linestyle=style_dict[stages_2[i]]['ls'])
        print(mockid_1)

        if 'RSD' in stages_2[i]:
            filename = utils.get_file_name(dirname,'gaussian-colore',N_side,pixel,compressed=True)
            h = fits.open(filename)
            lambdas_vel = lya*(1+h[4].data['Z'][:])
            i_col = np.where(h[1].data['MOCKID']==mockid_1)[0][0]
            print(h[1].data['MOCKID'][i_col])
            vel = h[3].data[i_col,:]
            axs[i].plot(lambdas_vel,vel*5000,label='vel*5000',color='grey',linestyle=style_dict[stages_2[i]]['ls'])
"""
