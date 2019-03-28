import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import glob

from pyacolore import plot_functions

default_location = '/global/homes/j/jfarr/Programs/picca/picca_analysis_000/'

if len(sys.argv) > 1:
    location = sys.argv[1]
else:
    location = default_location

linetypes = ['-','--',':','-:']
location_names = [r'$\beta=1.65$',r'$\beta=2.0$']

#Plotting decisions
afixs = ['free','fixed'] #'free' or 'fixed'
rmins = [20,40,60]
bb_colours = {'beta':'C0', 'bias':'C1', 'bias_eta':'C2', 'beta_BOSS':'C0', 'bias_BOSS':'C1', 'bias_eta_BOSS':'C2'}
a_colours = {'ap':'C0', 'at':'C1', 'ap_BOSS':'C0', 'at_BOSS':'C1', 'ap_combined':'C0', 'at_combined':'C1'}
show_plots = True
plot_combined = False

picca_runs = glob.glob(location+'/picca_*/')
if plot_combined:
    combined_runs = glob.glob(location+'/combined/')

#BOSS data from du Mas de Bourboux et al. 2017
z_BOSS = [2.4]
beta_BOSS = [1.650]
beta_BOSS_err = [0.081]
bias_BOSS = [-0.1337]
bias_BOSS_err = [0.00672]
bias_eta_BOSS = [-0.1337*1.650]
bias_eta_BOSS_err = [0.1337*1.650*np.sqrt((0.00672/0.1337)**2 + (0.081/1.650)**2)]
ap_BOSS = [1.069]
ap_BOSS_err = [0.029]
at_BOSS = [0.920]
at_BOSS_err = [0.034]

for afix in afixs:
    for rmin in rmins:
        all_bb_data = {}
        all_a_data = {}

        suffix = '_{}r_a{}'.format(rmin,afix)
        res_name = '/result'+suffix+'.h5'

        #Get the correlation objects
        corr_objects = plot_functions.get_correlation_objects(picca_runs,res_name=res_name)
        if plot_combined:
            combined_corr_object = plot_functions.get_correlation_objects(combined_runs,res_name=res_name)[0]

        bb_data = []
        a_data = []

        for corr_object in corr_objects:
            z_value = corr_object.zeff
            beta_LYA_value = corr_object.beta_LYA
            beta_LYA_error = corr_object.beta_LYA_err
            bias_LYA_value = corr_object.bias_LYA
            bias_LYA_error = corr_object.bias_LYA_err
            bias_LYA_eta_value = corr_object.bias_LYA_eta
            bias_LYA_eta_error = corr_object.bias_LYA_eta_err
            bb_data += [(z_value,beta_LYA_value,beta_LYA_error,bias_LYA_value,bias_LYA_error,bias_LYA_eta_value,bias_LYA_eta_error)]

            ap = corr_object.ap
            ap_err = corr_object.ap_err
            at = corr_object.at
            at_err = corr_object.at_err
            a_data += [(z_value,ap,ap_err,at,at_err)]
           
            chi2 = corr_object.fval
            ndat = corr_object.ndata
            npar = corr_object.npar
     
            print('|| {} || {:2.3f} +/- {:1.4f} || {:2.3f} +/- {:1.4f} || {:2.3f} +/- {:1.4f} || {:4.1f} / ({:4.0f} - {:1.0f})'.format(z_value,beta_LYA_value,beta_LYA_error,bias_LYA_value,bias_LYA_error,bias_LYA_eta_value,bias_LYA_eta_error,chi2,ndat,npar))
            #print('|| {} || {:2.3f} +/- {:1.4f} || {:2.3f} +/- {:1.4f} || {:4.1f} / ({:4.0f} - {:1.0f})'.format(z_value,ap,ap_err,at,at_err,chi2,ndat,npar))
 
        bb_dtype = [('z', 'd'), ('beta', 'd'), ('beta_err', 'd'), ('bias', 'd'), ('bias_err', 'd'), ('bias_eta', 'd'), ('bias_eta_err', 'd')]
        bb_data = np.array(bb_data,dtype=bb_dtype)
        bb_data = np.sort(bb_data,order=['z'])

        print(' ')
        a_dtype = [('z', 'd'), ('ap', 'd'), ('ap_err', 'd'), ('at', 'd'), ('at_err', 'd')]
        a_data = np.array(a_data,dtype=a_dtype)
        a_data = np.sort(a_data,order=['z'])

        if plot_combined:
            combined_z_value = combined_corr_object.zeff    
            combined_ap = combined_corr_object.ap
            combined_ap_err = combined_corr_object.ap_err
            combined_at = combined_corr_object.at
            combined_at_err = combined_corr_object.at_err

        #Plot bias, beta and bias_eta
        plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
        plt.errorbar(bb_data['z'],bb_data['beta'],yerr=bb_data['beta_err'],marker='o',label='beta(z)',color=bb_colours['beta'])
        plt.errorbar(bb_data['z'],bb_data['bias'],yerr=bb_data['bias_err'],marker='o',label='bias(z)',color=bb_colours['bias'])
        plt.errorbar(bb_data['z'],bb_data['bias_eta'],yerr=bb_data['bias_eta_err'],marker='o',label='bias_eta(z)',color=bb_colours['bias_eta'])
        plt.errorbar(z_BOSS,beta_BOSS,yerr=beta_BOSS_err,marker='x',label='beta BOSS',color=bb_colours['beta_BOSS'])
        plt.errorbar(z_BOSS,bias_BOSS,yerr=bias_BOSS_err,marker='x',label='bias BOSS',color=bb_colours['bias_BOSS'])
        plt.errorbar(z_BOSS,bias_eta_BOSS,yerr=bias_eta_BOSS_err,marker='x',label='bias_eta BOSS',color=bb_colours['bias_eta_BOSS'])        
        plt.legend()
        plt.grid()
        plt.xlabel('z')
        plt.title('Rmin = {} Mpc/h, alphas {}'.format(rmin,afix))
        plt.savefig(location+'/beta_biases'+suffix+'.pdf')
        if show_plots:
            plt.show()

        #Plot bias and bias_eta
        plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
        plt.errorbar(bb_data['z'],bb_data['bias'],yerr=bb_data['bias_err'],marker='o',label='bias(z)',color=bb_colours['bias'])
        plt.errorbar(bb_data['z'],bb_data['bias_eta'],yerr=bb_data['bias_eta_err'],marker='o',label='bias_eta(z)',color=bb_colours['bias_eta'])
        plt.errorbar(z_BOSS,bias_BOSS,yerr=bias_BOSS_err,marker='x',label='bias BOSS',color=bb_colours['bias_BOSS'])
        plt.errorbar(z_BOSS,bias_eta_BOSS,yerr=bias_eta_BOSS_err,marker='x',label='bias_eta BOSS',color=bb_colours['bias_eta_BOSS'])
        plt.legend()
        plt.grid()
        plt.xlabel('z')
        plt.title('Rmin = {} Mpc/h, alphas {}'.format(rmin,afix))
        plt.savefig(location+'/biases'+suffix+'.pdf')
        if show_plots:
            plt.show()
        
        #Plot alphas
        if afix == 'free':
            plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
            plt.errorbar(a_data['z'],a_data['ap'],yerr=a_data['ap_err'],marker='o',label='ap(z)')
            plt.errorbar(a_data['z'],a_data['at'],yerr=a_data['at_err'],marker='o',label='at(z)')
            plt.errorbar(z_BOSS,ap_BOSS,yerr=ap_BOSS_err,marker='x',label='ap BOSS',color=a_colours['ap_BOSS'])
            plt.errorbar(z_BOSS,at_BOSS,yerr=at_BOSS_err,marker='x',label='at BOSS',color=a_colours['at_BOSS'])
            if plot_combined:
                plt.errorbar(combined_z_value,combined_ap,yerr=combined_ap_err,marker='X',label='ap combined',color=a_colours['ap_combined'])
                plt.errorbar(combined_z_value,combined_at,yerr=combined_at_err,marker='X',label='at combined',color=a_colours['at_combined'])
            plt.legend()
            plt.grid()
            plt.xlabel('z')
            plt.title('Rmin = {} Mpc/h, alphas {}'.format(rmin,afix))
            plt.savefig(location+'/alphas'+suffix+'.pdf')
            if show_plots:
                plt.show()
