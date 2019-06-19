import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from astropy.io import fits
from picca import wedgize
import sys
from numpy import linalg
import mcfit
from lyacolore import correlation_model
import glob
import h5py
import warnings

class picca_correlation:
    def __init__(self,loc,p,cd,f):

        self.location = loc

        #picca parameters
        self.correl_type = p['correl_type']
        self.N_side = int(p['N_side'])
        self.N_pixels = int(p['N_pixels'])
        self.quantity_1 = p['quantities'][0]
        self.quantity_2 = p['quantities'][1]
        self.rpmin = float(p['rpmin'])
        self.rpmax = float(p['rpmax'])
        self.rtmin = float(p['rtmin'])
        self.rtmax = float(p['rtmax'])
        self.np = int(p['np'])
        self.nt = int(p['nt'])
        self.zmin = float(p['zmin'])
        self.zmax = float(p['zmax'])

        """
        #input skewer parameters
        self.sr = sr
        self.bm = bm
        """

        #correlation data in grids
        self.rp = cd['rp']
        self.rt = cd['rt']
        self.z = cd['z']
        self.xi = cd['xi']
        self.cov = cd['cov']
        self.nb = cd['nb']

        self.rp_grid = self.rp.reshape((self.np,self.nt))
        self.rt_grid = self.rt.reshape((self.np,self.nt))
        self.z_grid = self.z.reshape((self.np,self.nt))
        self.xi_grid = self.xi.reshape((self.np,self.nt))
        self.cov_grid = self.cov.reshape((self.np*self.nt,self.np*self.nt))
        self.nb_grid = self.nb.reshape((self.np,self.nt))

        """
        #cosmology data
        self.Z_cosmology = Z_cosmology
        self.D_cosmology = D_cosmology
        """

        if f:
            print(f.keys())
            #fit parameters
            self.zeff = f['zeff']
            self.fval = f['fval']
            self.ndata = f['ndata']
            self.npar = f['npar']

            self.growth_rate = f['growth_rate']['value']
            self.growth_rate_err = f['growth_rate']['error']

            if self.correl_type == 'cf' or self.correl_type == 'xcf':
                self.beta_LYA = f['beta_LYA']['value']
                self.beta_LYA_err = f['beta_LYA']['error']

                self.bias_LYA_eta = f['bias_eta_LYA']['value']
                self.bias_LYA_eta_err = f['bias_eta_LYA']['error']

                self.bias_LYA = self.bias_LYA_eta * self.growth_rate / self.beta_LYA
                self.bias_LYA_err = abs(self.bias_LYA * np.sqrt((self.bias_LYA_eta_err/self.bias_LYA_eta)**2 + (self.beta_LYA_err/self.beta_LYA)**2 + (self.growth_rate_err/self.growth_rate)**2))

            if self.correl_type == 'xcf' or self.correl_type == 'co':
                self.beta_QSO = f['beta_QSO']['value']
                self.beta_QSO_err = f['beta_QSO']['error']

                self.bias_QSO_eta = f['bias_eta_QSO']['value']
                self.bias_QSO_eta_err = f['bias_eta_QSO']['error']

                self.bias_QSO = self.bias_QSO_eta * self.growth_rate / self.beta_QSO
                self.bias_QSO_err = abs(self.bias_LYA * np.sqrt((self.bias_QSO_eta_err/self.bias_QSO_eta)**2 + (self.beta_QSO_err/self.beta_QSO)**2 + (self.growth_rate_err/self.growth_rate)**2))

            self.ap = f['ap']['value']
            self.ap_err = f['ap']['error']
            self.at = f['at']['value']
            self.at_err = f['at']['error']

            self.fit_xi_grid = f['xi_grid']

        return

    def get_biases_and_betas(self):

        if self.correl_type == 'cf':
            bias1 = self.bias_LYA
            bias2 = self.bias_LYA
            beta1 = self.beta_LYA
            beta2 = self.beta_LYA
        elif self.correl_type == 'xcf':
            bias1 = self.bias_LYA
            bias2 = self.bias_QSO
            beta1 = self.beta_LYA
            beta2 = self.beta_QSO
        elif self.correl_type == 'co':
            bias1 = self.bias_QSO
            bias2 = self.bias_QSO
            beta1 = self.beta_QSO
            beta2 = self.beta_QSO

        return bias1,bias2,beta1,beta2

    @classmethod
    def make_correlation_object(cls,location,file_name,res_name=None):

        #get parameters
        #replace with get_parameters_from_param_file when written
        parameters = get_parameters_from_param_file(location+'/parameters.txt')

        #get correlation data
        correlation_data = get_correlation_data(location+file_name)

        #get cosmology data
        #yet to do this

        if res_name:
            #get fit paramters
            fit_parameters = get_fit_from_result(location+res_name)
        else:
            fit_parameters = None

        return cls(location,parameters,correlation_data,fit_parameters)

    def plot_wedge(self,ax,mubin,plot_label,colour,r_power=2,nr=40,rmax=160.):

        mumin = mubin[0]
        mumax = mubin[1]

        b = wedgize.wedge(mumin=mumin,mumax=mumax,rtmax=self.rtmax,nrt=self.nt,
            rpmin=self.rpmin,rpmax=self.rpmax,nrp=self.np,nr=nr,rmax=rmax)

        """
        #REMOVE NON_DIAGONALS FROM COV MATRIX
        for i in range(self.cov_grid.shape[0]):
            for j in range(self.cov_grid.shape[1]):
                if i!=j:
                    self.cov_grid[i,j]=0
        """

        r, xi_wed, cov_wed = b.wedge(self.xi,self.cov)

        Nr = len(r)
        err_wed = np.zeros(Nr)
        for i in range(Nr):
            err_wed[i] = np.sqrt(cov_wed[i][i])

        cut = err_wed>0

        r = r[cut]
        xi_wed = xi_wed[cut]
        err_wed = err_wed[cut]

        ax.errorbar(r,(r**r_power) * xi_wed,yerr=(r**r_power) * err_wed,fmt='o',label=plot_label)

        return

    def plot_grid(self,plot_label,r_power,vmax=10**-4,xlabel='',ylabel='',label_fontsize=12,show_grid=True):

        fig = plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')

        im_grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(1,2),
                 axes_pad=0.15,
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="7%",
                 cbar_pad=0.15,
                 )

        #Add subplot containing the data.
        ax = im_grid[0]

        r_grid = np.sqrt(self.rp_grid**2 + self.rt_grid**2)

        to_show = self.xi_grid * (r_grid**r_power)

        #Mask the areas we don't want to plot
        if vmax:
            mask = to_show>vmax
        else:
            mask = np.zeros(to_show.shape)
        to_show = np.ma.masked_array(to_show,mask=mask)

        cmap = plt.cm.get_cmap('viridis', 20)

        im = ax.imshow(to_show,aspect='auto',cmap=cmap,origin='lower',
                    vmax=vmax,extent=[min(self.rt),max(self.rt),min(self.rp),
                    max(self.rp)])

        ax.set_title('Measured',fontsize=label_fontsize)
        ax.set_xlabel(xlabel,fontsize=label_fontsize)
        ax.set_ylabel(ylabel,fontsize=label_fontsize)
        if show_grid:
            ax.grid()

        #Add subplot showing the model.
        ax = im_grid[1]

        fit_grid = self.fit_xi_grid.reshape((self.np,self.nt))
        fit_to_show = fit_grid * (r_grid**r_power)

        #Mask the areas we don't want to plot.
        if vmax:
            mask = fit_to_show>vmax
        else:
            mask = np.zeros(fit_to_show.shape)
        fit_to_show = np.ma.masked_array(fit_to_show,mask=mask)
        im = ax.imshow(fit_to_show,aspect='auto',cmap=cmap,origin='lower',
                    vmax=vmax,extent=[min(self.rt),max(self.rt),min(self.rp),
                    max(self.rp)])

        ax.set_title('Fit',fontsize=label_fontsize)
        ax.set_xlabel(xlabel,fontsize=label_fontsize)
        ax.set_ylabel(ylabel,fontsize=label_fontsize)
        if show_grid:
            ax.grid()

        ax.cax.colorbar(im)
        ax.cax.toggle_label(True)

        return

    def plot_fit(self,ax,mubin,plot_label,colour,r_power=2,nr=40,rmax=160.):

        #Using data in picca fit
        mumin = mubin[0]
        mumax = mubin[1]

        b = wedgize.wedge(mumin=mumin,mumax=mumax,rtmax=self.rtmax,nrt=self.nt,
            rpmin=self.rpmin,rpmax=self.rpmax,nrp=self.np,nr=nr,rmax=rmax)

        """
        #REMOVE NON_DIAGONALS FROM COV MATRIX
        for i in range(self.cov_grid.shape[0]):
            for j in range(self.cov_grid.shape[1]):
                if i!=j:
                    self.cov_grid[i,j]=0
        """

        r, fit_xi_wed, _ = b.wedge(self.fit_xi_grid,self.cov)

        ax.plot(r,(r**r_power) * fit_xi_wed,c=colour,label=plot_label)

        return

    def plot_manual_model(self,ax,b1,b2,beta1,beta2,mubin,plot_label,colour,r_power=2,nr=40,rmax=160.):

        #Using data in picca fit
        mumin = mubin[0]
        mumax = mubin[1]

        b = wedgize.wedge(mumin=mumin,mumax=mumax,rtmax=self.rtmax,nrt=self.nt,
            rpmin=self.rpmin,rpmax=self.rpmax,nrp=self.np,nr=nr,rmax=rmax)

        #Using our model
        model = 'Slosar11'
        z_value = (self.zmin + self.zmax) / 2.
        r,xi = correlation_model.get_model_xi(model,self.quantity_1,self.quantity_2,b1,b2,beta1,beta2,z_value,mubin)
        indices = r<rmax
        r = r[indices]
        xi = xi[indices]

        ax.plot(r,(r**r_power) * xi,c=colour,label=plot_label)

        return

    def plot_vs_rt(self,np_bins,b1,b2,beta1,beta2,bin_list=None,r_power=0.,add_model=True):

        rp_vals_per_bin = self.np // np_bins

        if add_model:

            N_mult_res = 10
            np_model = self.np * N_mult_res
            nt_model = self.nt * N_mult_res
            rp_model = np.linspace(self.rpmin,self.rpmax,np_model)
            rt_model = np.linspace(self.rtmin,self.rtmax,nt_model)

            rp_model_grid = np.outer(rp_model,np.ones(nt_model))
            rt_model_grid = np.outer(np.ones(np_model),rt_model)

            z_value = (self.zmin + self.zmax) / 2.
            r_model_grid,xi_model_grid = correlation_model.get_model_xi_grid('Slosar11',self.quantity_1,self.quantity_2,b1,b2,beta1,beta2,z_value,rp_model_grid,rt_model_grid)

        for ind,i in enumerate(bin_list):
            col = 'C'+str(ind)

            rt = np.average(self.rt_grid[i*rp_vals_per_bin:(i+1)*rp_vals_per_bin,:], weights=self.nb_grid[i*rp_vals_per_bin:(i+1)*rp_vals_per_bin,:],axis=0)
            rp = np.average(self.rp_grid[i*rp_vals_per_bin:(i+1)*rp_vals_per_bin,:], weights=self.nb_grid[i*rp_vals_per_bin:(i+1)*rp_vals_per_bin,:],axis=0)

            xi = np.average(self.xi_grid[i*rp_vals_per_bin:(i+1)*rp_vals_per_bin,:], weights=self.nb_grid[i*rp_vals_per_bin:(i+1)*rp_vals_per_bin,:],axis=0)

            err_grid = np.sqrt(np.diag(self.cov_grid).reshape((self.np,self.nt)))

            xi_err = 1/np.sqrt(np.sum(1/(err_grid**2)[i*rp_vals_per_bin:(i+1)*rp_vals_per_bin,:],axis=0))

            r = np.sqrt(rp**2 + rt**2)
            xi_err_plot = xi_err*(r**r_power)
            xi_plot = xi*(r**r_power)

            drp = (self.rpmax - self.rpmin) / self.np

            plt.errorbar(rt,xi_plot,yerr=xi_err_plot,label=r'${:3.1f}<r_p<{:3.1f}$'.format(i*rp_vals_per_bin*drp,(i+1)*rp_vals_per_bin*drp),fmt='o',color=col)

            if add_model:

               rp_vals_per_bin_model = np_model // np_bins

               r_model = np.average(r_model_grid[i*rp_vals_per_bin_model:(i+1)*rp_vals_per_bin_model,:],axis=0)
               xi_model = np.average(xi_model_grid[i*rp_vals_per_bin_model:(i+1)*rp_vals_per_bin_model,:],axis=0)
               plt.plot(rt_model,xi_model*(r_model**r_power),c=col)
               #plt.plot(rt,xi_model*(r**r_power),c=col)

        return


def get_fit_from_result(filepath):

    ff = h5py.File(filepath,'r')
    fit = {}

    zeff = ff['best fit'].attrs['zeff']
    fit['zeff'] = zeff
    fval = ff['best fit'].attrs['fval']
    fit['fval'] = fval
    ndata = ff['best fit'].attrs['ndata']
    fit['ndata'] = ndata
    npar = ff['best fit'].attrs['npar']
    fit['npar'] = npar

    cosmo_pars = ["bias_eta_LYA","beta_LYA","bias_eta_QSO","beta_QSO","ap","at","growth_rate"]
    for par in cosmo_pars:
        if par in ff['best fit'].attrs:
            par_dict = {}
            value,error = ff['best fit'].attrs[par]
            par_dict['value'] = value
            par_dict['error'] = error
            fit[par] = par_dict
            print('{} = {} +/- {}'.format(par,value,error))

    # TODO: This won't work if it's not the lya auto correlation
    try:
        fit['xi_grid'] = ff['LYA(LYA)xLYA(LYA)/fit'][...]
    except KeyError:
        fit['xi_grid'] = ff['LYA(LYA)xQSO/fit'][...]

    ff.close()

    return fit

def get_correlation_data(filepath):

    #Use the cf_exp file.
    h = fits.open(filepath)
    correlation_data = {}
    correlation_data['rp'] = h[1].data['RP']
    correlation_data['rt'] = h[1].data['RT']
    correlation_data['z'] = h[1].data['Z']
    correlation_data['xi'] = h[1].data['DA']
    correlation_data['cov'] = h[1].data['CO']
    correlation_data['nb'] = h[1].data['NB']
    h.close()

    return correlation_data

def get_parameters_from_param_file(filepath):

    text_file = open(filepath, "r")
    lines = text_file.read().split('\n')
    params = {}
    for line in lines[:-1]:
        param = line.split(' = ')
        name = param[0]
        val = param[1]
        params[name] = val

    return params

def get_correlation_object(plot_info):

    location = plot_info['location']
    filename = plot_info['filename']
    res_name = 'result_{}r_a{}.h5'.format(str(int(plot_info['picca_fit_data']['rmin'])),plot_info['picca_fit_data']['afix'])
    corr_object = picca_correlation.make_correlation_object(location,filename,res_name=res_name)

    return corr_object

def make_colour_plots(corr_objects,r_power=0.,vmax=10**-4,save_plots=True,show_plots=True,suffix=''):

    for i,corr_object in enumerate(corr_objects):
        xlabel = r'$r_t\ /\ Mpc/h$'
        ylabel = r'$r_p\ /\ Mpc/h$'

        corr_object.plot_grid('',r_power,vmax=vmax,xlabel=xlabel,ylabel=ylabel)

        if save_plots:
            plt.savefig(corr_object.location+'/cf_colour'+suffix+'.pdf')

        if show_plots:
            plt.show()

    return

def plot_wedges(ax,plot_info):

    #Unpack the plot_info dictionary.
    corr_obj = plot_info['corr_object']

    for i,mubin in enumerate(plot_info['mu_bins']):
        #Plot the data.
        plot_label = r'${}<\mu<{}$'.format(mubin[0],mubin[1])
        colour = plot_info['mu_bin_colours'][i]
        corr_obj.plot_wedge(ax,mubin,plot_label,colour,**plot_info['plot_data'])

        #Add a model or fit.
        if plot_info['plot_picca_fit']:
            corr_obj.plot_fit(ax,mubin,'',colour,**plot_info['plot_data'])
        elif plot_info['plot_manual_fit']:
            b1,b2,beta1,beta2 = plot_info['manual_fit_data'].values()
            corr_obj.plot_fit(ax,mubin,b1,b2,beta1,beta2,'',colour,**plot_info['plot_data'])

    #Add axis labels.
    if plot_info['format']['xlabel']:
        ax.set_xlabel(r'$r\ [Mpc/h]$')
    if plot_info['format']['ylabel']:
        if plot_info['plot_data']['r_power'] > 1:
            ax.set_ylabel(r'$r^{} \xi (r)\ [(Mpc/h)^{}]$'.format(int(plot_info['plot_data']['r_power']),int(plot_info['plot_data']['r_power'])))
        elif plot_info['plot_data']['r_power'] == 1:
            ax.set_ylabel(r'$r\xi (r)\ [Mpc/h]$')
        elif plot_info['plot_data']['r_power'] == 0:
            ax.set_ylabel(r'$\xi (r)$')

    #Add a legend if desired.
    if plot_info['format']['legend']:
        ax.legend()

    return

def make_wedge_plots(corr_objects,mu_boundaries,plot_system,r_power,fit_type='picca',fit_data=None,nr=40,rmin=10.,rmax=160.,save_plots=True,show_plots=True,save_loc='.',suffix=''):

    colours = ['C0','C1','C2','C3','C4','C5','C6']

    mu_bins = bins_from_boundaries(mu_boundaries)

    if plot_system == 'plot_per_bin':

        colours = colours[:len(corr_objects)]
        for mu_bin in mu_bins:
            plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')

            for i,corr_object in enumerate(corr_objects):
                corr_object.plot_wedge(mu_bin,'',r_power,colours[i],nr=nr,rmax=rmax)
                if fit_type == 'picca':
                    corr_object.plot_fit(mu_bin,'',r_power,colours[i],nr=nr,rmax=rmax)
                elif fit_type == 'manual':
                    b1 = fit_data['b1']
                    b2 = fit_data['b2']
                    beta1 = fit_data['beta1']
                    beta2 = fit_data['beta2']
                    corr_object.plot_manual_model(b1,b2,beta1,beta2,mu_bin,'',r_power,colours[i],nr=nr,rmax=rmax)

            plt.title('{} < mu < {}'.format(mu_bin[0],mu_bin[1]))
            plt.legend(fontsize=12)
            plt.grid()

    elif plot_system == 'plot_per_file':

        colours = colours[:len(mu_bins)]
        for corr_object in corr_objects:
            plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')

            for i,mu_bin in enumerate(mu_bins):
                plot_label = r'{}<$\mu$<{}'.format(mu_bin[0],mu_bin[1])
                corr_object.plot_wedge(mu_bin,plot_label,r_power,colours[i])
                if fit_type == 'picca':
                    plot_label += ' (fit)'
                    corr_object.plot_fit(mu_bin,None,r_power,colours[i],nr=nr,rmax=rmax)

                    title_line_1 = r'{} {}{}; {} pix @ Nside {}; $z_{{eff}} = {:1.2f}$; $r_{{min}} = $??'.format(corr_object.correl_type,corr_object.quantity_1,corr_object.quantity_2,corr_object.N_pixels,16,corr_object.zeff)
                    title_line_2 = r'$b_\delta = {:1.3f}\pm{:1.3f}$; $\beta = {:1.3f}\pm{:1.3f}$; $\chi^2/(n_d-n_p) = {:5.1f}/({}-{})$; $\alpha_p = {:1.3f}\pm{:1.3f}$, $\alpha_t = {:1.3f}\pm{:1.3f}$'.format(corr_object.bias_LYA,corr_object.bias_LYA_err,corr_object.beta_LYA,corr_object.beta_LYA_err,corr_object.fval,corr_object.ndata,corr_object.npar,corr_object.ap,corr_object.ap_err,corr_object.at,corr_object.at_err)

                elif fit_type == 'manual':
                    plot_label += ' (model)'
                    b1 = fit_data['b1']
                    b2 = fit_data['b2']
                    beta1 = fit_data['beta1']
                    beta2 = fit_data['beta2']
                    corr_object.plot_manual_model(b1,b2,beta1,beta2,mu_bin,plot_label,r_power,colours[i],nr=nr,rmax=rmax)

                    title_line_1 = r'{} {}{}; {} pix @ Nside {}; ${} < z < {}$; $r_{{min}} = $??'.format(corr_object.correl_type,corr_object.quantity_1,corr_object.quantity_2,corr_object.N_pixels,16,corr_object.zmin,corr_object.zmax)
                    title_line_2 = r'Model is Legendre decomposition using: $b_\delta$s$ = {:1.3f},{:1.3f}$; $\beta$s$ = {:1.3f},{:1.3f}$'.format(b1,b2,beta1,beta2)

            title = title_line_1 + '\n' + title_line_2

            # TODO: import the LyaCoLoRe nside to use here, use corr_object.correl_type to determine which biases are presented
            plt.title(title)
            plt.legend(fontsize=12)
            plt.grid()
            plt.xlabel(r'$r\ /\ Mpc/h$',fontsize=12)
            plt.ylabel(r'$r^2 \xi (r)$',fontsize=12)

            if save_plots:
                plt.savefig(corr_object.location+'/cf'+suffix+'_'+fit_type+'fit.pdf')

    if show_plots:
        plt.show()

    return

def bins_from_boundaries(boundaries):

    bins = list(zip(boundaries[:-1],boundaries[1:]))

    return bins

def make_plot_vs_rt(corr_objects,np_bins,fit_data,bin_list=None,r_power=0.,show_plots=True,save_plots=True):

    for corr_object in corr_objects:
        plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
        b1 = fit_data['b1']
        b2 = fit_data['b2']
        beta1 = fit_data['beta1']
        beta2 = fit_data['beta2']
        corr_object.plot_vs_rt(np_bins,b1,b2,beta1,beta2,bin_list=bin_list,r_power=r_power)
        plt.grid()
        plt.legend()
        plt.xlabel(r'$r_t$')
        plt.ylabel(r'$\xi r^{}$'.format(int(r_power)))

        title_line_1 = r'{} {}{}; {} pix @ Nside {}; ${} < z < {}$; $r_{{min}} = $??'.format(corr_object.correl_type,corr_object.quantity_1,corr_object.quantity_2,corr_object.N_pixels,16,corr_object.zmin,corr_object.zmax)
        title_line_2 = r'Model is Legendre decomposition using: $b_\delta$s$ = {:1.3f},{:1.3f}$; $\beta$s$ = {:1.3f},{:1.3f}$'.format(b1,b2,beta1,beta2)
        title = title_line_1 + '\n' + title_line_2
        plt.title(title)

        if save_plots:
            plt.savefig(corr_object.location+'/cf_xi_vs_rt.pdf')
        if show_plots:
            plt.show()

    return


###############################################################################




def get_scales(filepaths):

    scales = []
    for filepath in filepaths:
        scale = 1.0
        parameters = get_parameters_from_filename(filepath)
        z = (parameters['zmax']+parameters['zmin'])/2
        for quantity in parameters['quantities']:
            scale *= correlation_model.get_bias(z,quantity)
            scale *= correlation_model.get_growth_factor_scaling(z,quantity)
        scales += [scale]

    return scales

def get_parameters_from_filename(file_path):

    default_cf_parameters = {'correlation': 'cf', 'nside': 16, 'sr': 2.0, 'rpmax': 200.0, 'rpmin': 0.0, 'rtmax': 200.0, 'np': 50, 'nt': 50, 'nr': 25, 'rmax': 200.0, 'zmin': 0.0, 'zmax':4.0, 'quantities': 'GG', 'bm': '3'}
    default_xcf_parameters = {'correlation': 'xcf', 'nside': 16, 'sr': 2.0, 'rpmax': 200.0, 'rpmin': -200.0, 'rtmax': 200.0, 'np': 100, 'nt': 50, 'nr': 25, 'rmax': 200.0, 'zmin': 0.0, 'zmax': 4.0, 'quantities': 'Gq', 'bm': '3'}
    parameter_dtypes = {'correlation': 'str', 'nside': 'int', 'sr': 'float', 'rpmax': 'float', 'rpmin': 'float', 'rtmax': 'float', 'np': 'int', 'nt': 'int', 'nr': 'int', 'rmax': 'float', 'zmin': 'float', 'zmax': 'float', 'quantities': 'str', 'bm': 'int'}


    #default_cf_parameters = {'correlation': 'cf', 'nside': 16, 'sr': 2.0, 'rpmax': 160.0, 'rpmin': 0.0, 'rtmax': 160.0, 'np': 40, 'nt': 40, 'nr': 20, 'rmax': 160.0, 'zmin': 0.0, 'zmax':4.0, 'quantities': 'GG', 'bm': '3'}


    last_slash = file_path[::-1].find('/')
    filename = file_path[len(file_path)-last_slash:]

    if filename[:3] == 'xcf':
        parameters = default_xcf_parameters
    elif filename[:2] == 'cf':
        parameters = default_cf_parameters

    for parameter in parameters.keys():
        if parameter in filename:
            i_pvs = int(filename.find('_'+parameter)+len(parameter)+1)
            i_pve = int(i_pvs + filename[i_pvs:].find('_'))
            if i_pve < i_pvs:
                i_pve = int(i_pvs + filename[i_pvs:].find('.fits'))
            d_type = type(parameter)

            if parameter_dtypes[parameter] == 'float':
                parameters[parameter] = float(filename[i_pvs:i_pve])
            elif parameter_dtypes[parameter] == 'int':
                parameters[parameter] = int(filename[i_pvs:i_pve])
            elif parameter_dtypes[parameter] == 'str':
                parameters[parameter] = str(filename[i_pvs:i_pve])

    return parameters

def get_plot_data(mumin,mumax,file_path):

    #Read data from the filename to know how to wedgize etc
    #filename format: ${CORRELTYPE}_${NPIXELS}_${BRANCH}_${OPTION}_${RSDOPTION}_rpmin${RPMIN}_rpmax${RPMAX}_rtmax${RTMAX}_np${NP}_nt${NT}_zmin${ZQSOMIN}_zmax${ZQSOMAX}.fits.gz
    parameters = get_parameters_from_filename(file_path)

    #file_parameters = {**file_parameters,**{filename:parameters}}

    h = fits.open(file_path)
    data = h[1].data

    b = wedgize.wedge(mumin=mumin,mumax=mumax,rtmax=parameters['rtmax'],
        nrt=int(parameters['nt']),rpmin=parameters['rpmin'],
        rpmax=parameters['rpmax'],nrp=int(parameters['np']),
        nr=int(parameters['nr']),rmax=parameters['rmax'])

    rp = data['RP'][:]
    rt = data['RT'][:]
    z = data['Z'][:]
    xi_grid = data['DA'][:]
    cov_grid = data['CO'][:]


    #REMOVE NON_DIAGONALS FROM COV MATRIX
    for i in range(cov_grid.shape[0]):
        for j in range(cov_grid.shape[1]):
            if i!=j:
                cov_grid[i,j]=0

    print(xi_grid.shape,cov_grid.shape)

    r, xi_wed, cov_wed = b.wedge(xi_grid,cov_grid)

    Nr = len(r)
    err_wed = np.zeros(Nr)
    for i in range(Nr):
        err_wed[i] = np.sqrt(cov_wed[i][i])

    cut = err_wed>0

    return r, xi_wed, err_wed, cut, parameters

def get_CAMB_xi(filepath,scale_CAMB,CAMB_sr):

    data = np.loadtxt(filepath)
    xi = data[:,1]
    r = data[:,0]
    xi = xi*scale_CAMB
    plot_label='CAMB, scaling: {:2.2f}, sr: {}'.format(scale_CAMB,CAMB_sr)

    return r,xi,plot_label

## TODO: Do I need this?
def add_CAMB_Pk(location,filename,z,mu,RSDOPTION='NO_RSD',CAMB_input='xi'):
    data = np.loadtxt(location+filename)
    k = data[:,0]
    Pk = data[:,1]
    return

def plot_xir2(mubin,filename,plot_label):

    mumin = mubin[0]
    mumax = mubin[1]

    r,xi_wed,err_wed,cut,parameters = get_plot_data(mumin,mumax,filename)

    plt.errorbar(r[cut],xi_wed[cut]*(r[cut]**2),yerr=err_wed[cut]*(r[cut]**2),fmt='o',label=plot_label)

    return parameters

def plot_per_bin(mubins,filenames,add_CAMB,plot_label_parameters,CAMB_sr=None,scale_CAMB=None,CAMB_location=None,CAMB_filename_format=None):

    for mubin in mubins:

        mumin = mubin[0]
        mumax = mubin[1]

        plt.figure()
        plt.axhline(y=0,color='gray',ls=':')
        plt.xlabel('r [Mpc/h]')
        plt.ylabel('r^2 xi(r)')
        plt.grid(True, which='both')
        plt.title('correlation function, {} < mu < {}'.format(mumin,mumax))

        for filename in filenames:
            parameters = get_parameters_from_filename(filename)
            plot_label = ''
            for parameter in plot_label_parameters:
                plot_label += '{}: {}, '.format(parameter,parameters[parameter])
                plot_label = plot_label[:-2]
            plot_xir2(mubin,filename,plot_label)

        if add_CAMB == True:
            for i,sr in enumerate(CAMB_sr):
                r_CAMB,xi_CAMB,plot_label_CAMB = get_CAMB_xi(CAMB_location+CAMB_filename_format.format(sr),scale_CAMB[i],sr)
                plt.plot(r_CAMB,xi_CAMB*(r_CAMB**2),label=plot_label_CAMB)

        #save figure
        plt.legend(loc=3)
        plt.savefig('xi_wedge_{}_{}.pdf'.format(mumin,mumax))

def plot_per_file(mubins,filenames,add_CAMB,plot_label_parameters,CAMB_sr=None,scale_CAMB=None,CAMB_location=None,CAMB_filename_format=None):

    for filename in filenames:

        plt.figure()
        plt.axhline(y=0,color='gray',ls=':')
        plt.xlabel('r [Mpc/h]')
        plt.ylabel(r'$r^2 \xi(r)$')
        plt.grid(True, which='both')

        plot_title = ''
        parameters = get_parameters_from_filename(filename)
        for parameter in plot_label_parameters:
            plot_title += '{}: {}, '.format(parameter,parameters[parameter])
            plot_title = plot_title[:-2]
        plt.title('correlation function, {}'.format(plot_title))

        for mubin in mubins:
            mumin = mubin[0]
            mumax = mubin[1]
            plot_label = '{} < mu < {}'.format(mumin,mumax)
            plot_xir2(mubin,filename,plot_label)

        if add_CAMB == True:
            for i,sr in enumerate(CAMB_sr):
                r_CAMB,xi_CAMB,plot_label_CAMB = get_CAMB_xi(CAMB_location+CAMB_filename_format.format(sr),scale_CAMB[i],sr)
                plt.plot(r_CAMB,xi_CAMB*(r_CAMB**2),label=plot_label_CAMB)

        #save figure
        plt.legend(loc=3)
        plt.savefig('xi_wedges_{}.pdf'.format(plot_title))

def visual_fit(filename,b_values,beta_values,model,data_parameters,z,compute_b_beta=False):

    mubin_boundaries = [0.0,1.0]
    mubins = []
    for i in range(len(mubin_boundaries)-1):
        mubins += [(mubin_boundaries[i],mubin_boundaries[i+1])]
    #mubins = [(0.0,0.33),(0.33,0.67),(0.67,1.0)]
    mubins = [(0.0,0.5),(0.5,0.8),(0.8,1.0)]
    #find a more sophisticated way to do this
    colours = ['r',(0.5,0.5,0.5),'b']

    #mubins = [(0.0,1.0)]
    #colours = ['r']

    N_bins = len(mubins)

    plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')

    for i,bin in enumerate(mubins):
        mumin = bin[0]
        mumax = bin[1]

        #Wedgize the data
        r,xi_wed,err_wed,cut,_ = get_plot_data(mumin,mumax,filename)

        data_label = 'data, {}<$\mu$<{}'.format(mumin,mumax)
        plt.errorbar(r[cut],xi_wed[cut]*(r[cut]**2),yerr=err_wed[cut]*(r[cut]**2),fmt='o',label=data_label,c=colours[i])

    quantity1 = data_parameters['quantities'][0]
    quantity2 = data_parameters['quantities'][1]

    #TODO: implement compute_b_beta. Need to have different flags for each quantity?
    for b1 in b_values[quantity1]:
        for beta1 in beta_values[quantity1]:
            for b2 in b_values[quantity2]:
                for beta2 in beta_values[quantity2]:

                    r_model,xi_model_values = correlation_model.get_model_xi(model,[b1,b2],[beta1,beta2],data_parameters,z,mubins,b_from_z=False)

                    for i,key in enumerate(xi_model_values.keys()):
                        #model_label = 'b_{}={}, beta_{}={}, b_{}={}, beta_{}={}, mu={}'.format(quantity1,b1,quantity1,beta1,quantity2,b2,quantity2,beta2,key)
                        model_label = 'model, $\mu$={}'.format(key)
                        plt.plot(r_model,xi_model_values[key]*(r_model**2),label=model_label,c=colours[i])
                        #plt.plot(r[cut],(xi_wed[cut])/(np.interp(r[cut],r_model,xi_model_values[key])),label='(RATIO): '+model_label)

    plt.axhline(y=0,color='gray',ls=':')
    plt.xlabel(r'$r / Mpc h^{-1}$')
    plt.ylabel(r'$r^2 \xi(r)$')
    plt.grid()
    plt.legend()
    plt.xlim(0,200)
    #Make this more sensible filename
    plt.savefig('fit_picca_visual.pdf')
    plt.show()

    return
