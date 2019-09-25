import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from astropy.io import fits
from picca import wedgize
import sys
from numpy import linalg
import mcfit
from lyacolore import correlation_model, utils
import glob
import h5py
import warnings

class picca_correlation:
    def __init__(self,loc,ci,cd,f):

        self.location = loc

        #correlation information
        self.correl_type = ci['correl_type']
        try:
            self.quantity_1 = ci['quantities'][0]
            self.quantity_2 = ci['quantities'][1]
            for q in [self.quantity_1,self.quantity_2]:
                if q == 'q':
                    q == 'QSO'
                elif q == 'F':
                    q = 'LYA'
        except:
            self.quantity_1 = ci['quantity 1']
            self.quantity_2 = ci['quantity 2']

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

        self.rpmin = cd['rpmin']
        self.rpmax = cd['rpmax']
        self.rtmin = cd['rtmin']
        self.rtmax = cd['rtmax']
        self.np = cd['np']
        self.nt = cd['nt']

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
            self.fit = f
            self.zeff = f['zeff']
        else:
            self.zeff = utils.get_zeff(self.z,self.rp,self.rt,self.nb,rmin=80.,rmax=120.)

        print('zeff:',self.zeff)
        return

    @classmethod
    def make_correlation_object(cls,location,file_name,res_name=None,res_location=None,corr_name='LYA(LYA)xLYA(LYA)'):

        if res_location is None:
            res_location = location

        try:
            parameters = get_parameters_from_param_file(location+'/parameters.txt')
        except:
            parameters = get_parameters_from_param_file(location+'/corr_info.txt')

        #get correlation data
        correlation_data = get_correlation_data(location+file_name)

        #get cosmology data
        #yet to do this

        if res_name:
            #get fit paramters
            fit_parameters = get_fit_from_result(res_location,res_name,corr_name,parameters['correl_type'])
        else:
            fit_parameters = None

        return cls(location,parameters,correlation_data,fit_parameters)

    def plot_wedge(self,ax,mubin,plot_label,colour,r_power=2,nr=40,rmax_plot=200.,abs_mu=False):

        mumin = mubin[0]
        mumax = mubin[1]

        #Create the wedge, and wedgise the correlation.
        b = wedgize.wedge(mumin=mumin,mumax=mumax,rtmax=self.rtmax,nrt=self.nt,
            rpmin=self.rpmin,rpmax=self.rpmax,nrp=self.np,nr=nr,rmax=rmax_plot,
            absoluteMu=abs_mu)
        r, xi_wed, cov_wed = b.wedge(self.xi,self.cov)

        print(mubin)
        print(r[:4])
        print(xi_wed[:4])
        print(' ')

        #Get the errors.
        Nr = len(r)
        err_wed = np.zeros(Nr)
        for i in range(Nr):
            err_wed[i] = np.sqrt(cov_wed[i][i])

        #Define the variables to plot, and plot them.
        cut = err_wed>0
        r = r[cut]
        xi_wed = xi_wed[cut]
        err_wed = err_wed[cut]
        ar = ax.errorbar(r,(r**r_power) * xi_wed,yerr=(r**r_power) * err_wed,fmt='o',label=plot_label)

        return ar, plot_label

    def plot_wedge_fit(self,ax,mubin,plot_label,colour,r_power=2,nr=40,rmax_plot=200.,rmin_fit=None,rmax_fit=None,abs_mu=False):

        #Using data in picca fit
        mumin = mubin[0]
        mumax = mubin[1]

        b = wedgize.wedge(mumin=mumin,mumax=mumax,rtmax=self.rtmax,nrt=self.nt,
            rpmin=self.rpmin,rpmax=self.rpmax,nrp=self.np,nr=nr,rmax=rmax_plot,
            absoluteMu=abs_mu)

        """
        #REMOVE NON_DIAGONALS FROM COV MATRIX
        for i in range(self.cov_grid.shape[0]):
            for j in range(self.cov_grid.shape[1]):
                if i!=j:
                    self.cov_grid[i,j]=0
        """

        r, fit_xi_wed, _ = b.wedge(self.fit['xi_grid'],self.cov)

        if rmin_fit == None:
            rmin_fit = 0.
        if rmax_fit == None:
            rmax_fit = rmax_plot

        #Plot the fit up to rmax_plot in a dashed line.
        ax.plot(r,(r**r_power) * fit_xi_wed,c=colour,linestyle='--')

        #Plot the fit in the fitted region in a solid line.
        fit_bins = (r>rmin_fit) * (r<rmax_fit)
        ax.plot(r[fit_bins],(r[fit_bins]**r_power) * fit_xi_wed[fit_bins],c=colour,label=plot_label)

        return

    def plot_wedge_manual_model(self,ax,b1,b2,beta1,beta2,mubin,plot_label,colour,r_power=2,nr=40,rmax_plot=200.,abs_mu=False):

        #Using data in picca fit
        mumin = mubin[0]
        mumax = mubin[1]

        b = wedgize.wedge(mumin=mumin,mumax=mumax,rtmax=self.rtmax,nrt=self.nt,
            rpmin=self.rpmin,rpmax=self.rpmax,nrp=self.np,nr=nr,rmax=rmax_plot,
            absoluteMu=abs_mu)

        #Using our model
        model = 'Slosar11'
        r,xi = correlation_model.get_model_xi(model,self.quantity_1,self.quantity_2,b1,b2,beta1,beta2,self.zeff,mubin)
        indices = r<rmax_plot
        r = r[indices]
        xi = xi[indices]

        ax.plot(r,(r**r_power) * xi,c=colour,label=plot_label,linestyle=':')

        return

    def plot_rp_bin_vs_rt(self,ax,rpbin,plot_label,colour,r_power=2,nr=40,rmax_plot=200.):

        rpmin = rpbin[0]
        rpmax = rpbin[1]

        #Get the weights of how to group the bins in rp.
        rp_edges = np.linspace(self.rpmin,self.rpmax,self.np+1)
        rp_bin_widths = rp_edges[1:] - rp_edges[:-1]
        weights = np.maximum(0.,np.minimum(rp_bin_widths,np.minimum(rpmax-rp_edges[:-1],rp_edges[1:]-rpmin))) / rp_bin_widths
        weights_grid = np.outer(weights,np.ones(self.nt)).astype('int')

        rt = np.average(self.rt_grid,weights=weights_grid,axis=0)
        rp = np.average(self.rp_grid,weights=weights_grid,axis=0)
        xi = np.average(self.xi_grid,weights=weights_grid,axis=0)
        err_grid = np.sqrt(np.diag(self.cov_grid)).reshape((self.np,self.nt))

        # TODO: Not sure about this. Need to use Nb for when the number of bins is greater than 1.
        N_non_empty_rp_bins = np.sum(np.sum(weights_grid,axis=1)>0)
        if N_non_empty_rp_bins == 1:
            xi_err = err_grid[weights_grid>0]
        elif N_non_empty_rp_bins > 1:
            xi_err = 1 / np.sqrt(np.sum(1/((err_grid[weights_grid>0]*weights_grid[weights_grid>0])**2).reshape((N_non_empty_rp_bins,self.nt)),axis=0))

        #Define the variables to plot, and plot them.
        r = np.sqrt(rp**2 + rt**2)
        xi_err_plot = xi_err*(r**r_power)
        xi_plot = xi*(r**r_power)
        ax.errorbar(rt,xi_plot,yerr=xi_err_plot,label=r'${:3.1f}<r_p<{:3.1f}$'.format(rpmin,rpmax),fmt='o',color=colour)

        return

    def plot_rp_bin_vs_rt_fit(self,ax,rtbin,plot_label,colour,r_power=2,nr=40,rmax_plot=200.,rmin_fit=None,rmax_fit=None):

        #Using data in picca fit
        rtmin = rtbin[0]
        rtmax = rtbin[1]

        #Get the weights of how to group the bins in rp.
        rp_edges = np.linspace(self.rpmin,self.rpmax,self.np+1)
        rp_bin_widths = rp_edges[1:] - rp_edges[:-1]
        weights = np.maximum(0.,np.minimum(rp_bin_widths,np.minimum(rpmax-rp_edges[:-1],rp_edges[1:]-rpmin))) / rp_bin_widths
        weights_grid = np.outer(weights,np.ones(self.nt)).astype('int')

        rt = np.average(self.rt_grid,weights=weights_grid,axis=0)
        rp = np.average(self.rp_grid,weights=weights_grid,axis=0)
        xi_grid = self.fit['xi_grid'].reshape((self.np,self.nt))
        xi_model = np.average(xi_grid,weights=weights_grid,axis=0)

        #Plot the model.
        r = np.sqrt(rp**2 + rt**2)
        xi_plot = xi_model*(r**r_power)
        ax.plot(rp_model,xi_plot,c=colour)

        return

    def plot_rp_bin_vs_rt_manual_model(self,ax,b1,b2,beta1,beta2,rpbin,plot_label,colour,r_power=2,nr=40,rmax_plot=200.):

        #Using data in picca fit
        rpmin = rpbin[0]
        rpmax = rpbin[1]

        #Set up a finer rp-rt grid to plot our model with.
        model = 'Slosar11'
        N_mult_res = 10
        np_model = self.np * N_mult_res
        nt_model = self.nt * N_mult_res
        rp_model_edges = np.linspace(self.rpmin,self.rpmax,np_model+1)
        rt_model_edges = np.linspace(self.rtmin,self.rtmax,nt_model+1)
        rp_model = utils.get_centres(rp_model_edges)
        rt_model = utils.get_centres(rt_model_edges)
        rp_model_grid = np.outer(rp_model,np.ones(nt_model))
        rt_model_grid = np.outer(np.ones(np_model),rt_model)
        r_model_grid,xi_model_grid = correlation_model.get_model_xi_grid(model,self.quantity_1,self.quantity_2,b1,b2,beta1,beta2,self.zeff,rp_model_grid,rt_model_grid,sr=0.0)

        #Calculate the weights on this finer grid.
        rp_bin_widths = rp_model_edges[1:] - rp_model_edges[:-1]
        weights = np.maximum(0.,np.minimum(rp_bin_widths,np.minimum(self.rpmax-rp_model_edges[:-1],rp_model_edges[1:]-self.rpmin))) / rp_bin_widths
        weights_grid = np.outer(weights,np.ones(nt_model)).astype('int')

        print(r_model_grid.shape,xi_model_grid.shape)
        #Plot the model.
        r_model = np.average(r_model_grid,weights=weights_grid,axis=0)
        xi_model = np.average(xi_model_grid,weights=weights_grid,axis=0)
        ax.plot(rt_model,xi_model*(r_model**r_power),c=colour)

        return

    def plot_rt_bin_vs_rp(self,ax,rtbin,plot_label,colour,r_power=2,nr=40,rmax=160.,rmax_plot=200.):

        rtmin = rtbin[0]
        rtmax = rtbin[1]

        #Get the weights of how to group the bins in rt.
        rt_edges = np.linspace(self.rtmin,self.rtmax,self.nt+1)
        rt_bin_widths = rt_edges[1:] - rt_edges[:-1]
        weights = np.maximum(0.,np.minimum(rt_bin_widths,np.minimum(rtmax-rt_edges[:-1],rt_edges[1:]-rtmin))) / rt_bin_widths
        weights_grid = np.outer(np.ones(self.np),weights).astype('int')

        rt = np.average(self.rt_grid,weights=weights_grid,axis=1)
        rp = np.average(self.rp_grid,weights=weights_grid,axis=1)
        xi = np.average(self.xi_grid,weights=weights_grid,axis=1)
        err_grid = np.sqrt(np.diag(self.cov_grid)).reshape((self.np,self.nt))

        # TODO: Not sure about this. Need to use Nb for when the number of bins is greater than 1.
        N_non_empty_rt_bins = np.sum(np.sum(weights_grid,axis=0)>0)
        if N_non_empty_rt_bins == 1:
            xi_err = err_grid[weights_grid>0]
        elif N_non_empty_rt_bins > 1:
            xi_err = 1 / np.sqrt(np.sum(1/((err_grid[weights_grid>0]*weights_grid[weights_grid>0])**2).reshape((self.np,N_non_empty_rt_bins)),axis=1))

        #Define the variables to plot, and plot them.
        r = np.sqrt(rp**2 + rt**2)
        xi_err_plot = xi_err*(r**r_power)
        xi_plot = xi*(r**r_power)
        ax.errorbar(rp,xi_plot,yerr=xi_err_plot,label=r'${:3.1f}<r_\perp<{:3.1f}$'.format(rtmin,rtmax),fmt='o',color=colour)

        #Hack to plot residual
        #xi_model_grid = self.fit['xi_grid'].reshape((self.np,self.nt))
        #xi_model = np.average(xi_model_grid,weights=weights_grid,axis=1)
        #xi_model_plot = xi_model*(r**r_power)
        #ax.errorbar(rp,(xi_plot-xi_model_plot)/xi_model_plot,yerr=xi_err_plot,label=r'${:3.1f}<r_\perp<{:3.1f}$'.format(rtmin,rtmax),fmt='o',color=colour)


        return

    def plot_rt_bin_vs_rp_fit(self,ax,rtbin,plot_label,colour,r_power=2,nr=40,rmax_plot=200.,rmin_fit=None,rmax_fit=None):

        rtmin = rtbin[0]
        rtmax = rtbin[1]

        #Get the weights of how to group the bins in rt.
        rt_edges = np.linspace(self.rtmin,self.rtmax,self.nt+1)
        rt_bin_widths = rt_edges[1:] - rt_edges[:-1]
        weights = np.maximum(0.,np.minimum(rt_bin_widths,np.minimum(rtmax-rt_edges[:-1],rt_edges[1:]-rtmin))) / rt_bin_widths
        weights_grid = np.outer(np.ones(self.np),weights).astype('int')

        rt = np.average(self.rt_grid,weights=weights_grid,axis=1)
        rp = np.average(self.rp_grid,weights=weights_grid,axis=1)
        xi_grid = self.fit['xi_grid'].reshape((self.np,self.nt))
        xi = np.average(xi_grid,weights=weights_grid,axis=1)

        #Define the variables to plot, and plot them.
        r = np.sqrt(rp**2 + rt**2)
        xi_plot = xi*(r**r_power)
        ax.plot(rp,xi_plot,c=colour)

        return

    def plot_rt_bin_vs_rp_manual_model(self,ax,b1,b2,beta1,beta2,rtbin,plot_label,colour,r_power=2,nr=40,rmax_plot=200.):

        #Using data in picca fit
        rtmin = rtbin[0]
        rtmax = rtbin[1]

        #Set up a finer rp-rt grid to plot our model with.
        model = 'Slosar11'
        N_mult_res = 1
        np_model = self.np * N_mult_res
        nt_model = self.nt * N_mult_res
        rp_model_edges = np.linspace(self.rpmin,self.rpmax,np_model+1)
        rt_model_edges = np.linspace(self.rtmin,self.rtmax,nt_model+1)
        rp_model = utils.get_centres(rp_model_edges)
        rt_model = utils.get_centres(rt_model_edges)
        rp_model_grid = np.outer(rp_model,np.ones(nt_model))
        rt_model_grid = np.outer(np.ones(np_model),rt_model)
        r_model_grid,xi_model_grid = correlation_model.get_model_xi_grid(model,self.quantity_1,self.quantity_2,b1,b2,beta1,beta2,self.zeff,rp_model_grid,rt_model_grid,sr=0.0)

        #Calculate the weights on this finer grid.
        rt_bin_widths = rt_model_edges[1:] - rt_model_edges[:-1]
        weights = np.maximum(0.,np.minimum(rt_bin_widths,np.minimum(rtmax-rt_model_edges[:-1],rt_model_edges[1:]-rtmin))) / rt_bin_widths
        weights_grid = np.outer(np.ones(np_model),weights).astype('int')

        #Plot the model.
        r_model = np.average(r_model_grid,weights=weights_grid,axis=1)
        xi_model = np.average(xi_model_grid,weights=weights_grid,axis=1)
        ax.plot(rp_model,xi_model*(r_model**r_power),c=colour)

        return

    def plot_grid(self,ax,mubin,plot_label,colour,r_power=2,nr=40,rmax_plot=160.):
    #def plot_grid(self,plot_label,r_power,vmax=10**-4,xlabel='',ylabel='',label_fontsize=12,show_grid=True):

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

        fit_grid = self.fit['xi_grid'].reshape((self.np,self.nt))
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


def get_fit_from_result(location,result_name,corr_name,corr_type):

    result_filepath = location + result_name
    ff = h5py.File(result_filepath,'r')
    fit = {}
    for attr in ff['best fit'].attrs:
        p = ff['best fit'].attrs[attr]
        if isinstance(p,np.ndarray):
            par_dict = {}
            par_dict['value'] = p[0]
            par_dict['error'] = p[1]
            fit[attr] = par_dict
        else:
            fit[attr] = p
        #print(attr,par_dict)

    print('chi2: {:4.1f}/({:d}-{:d})'.format(fit['fval'],fit['ndata'],fit['npar']))

    if corr_name:
        fit['xi_grid'] = ff['{}/fit'.format(corr_name)][...]
    else:
        if corr_type == 'cf':
            fit['xi_grid'] = ff['LYA(LYA)xLYA(LYA)/fit'][...]
        elif corr_type == 'xcf':
            fit['xi_grid'] = ff['LYA(LYA)xQSO/fit'][...]
        elif corr_type == 'co':
            fit['xi_grid'] = ff['QSOxQSO/fit'][...]

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

    correlation_data['rpmin'] = h[1].header['RPMIN']
    correlation_data['rpmax'] = h[1].header['RPMAX']
    try:
        correlation_data['rtmin'] = h[1].header['RTMIN']
    except:
        correlation_data['rtmin'] = 0.0
    correlation_data['rtmax'] = h[1].header['RTMAX']
    correlation_data['np'] = h[1].header['NP']
    correlation_data['nt'] = h[1].header['NT']

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
    if plot_info['plot_picca_fit']:
        try:
            res_name = plot_info['result_name']
        except KeyError:
            res_name = 'result_{}r_a{}.h5'.format(str(int(plot_info['picca_fit_data']['rmin'])),plot_info['picca_fit_data']['afix'])
        try:
            res_location = plot_info['result_location']
        except KeyError:
            res_location = location
    else:
        res_name = None
        res_location = None

    try:
        corr_name = plot_info['corr_name']
    except KeyError:
        corr_name = None

    corr_object = picca_correlation.make_correlation_object(location,filename,res_name=res_name,res_location=res_location,corr_name=corr_name)

    return corr_object

def make_colour_plots(ax,plot_info,vmax=10**-4):

    #Unpack the plot_info dictionary.
    corr_obj = plot_info['corr_object']

    for i,corr_object in enumerate(corr_objects):
        xlabel = r'$r_t\ /\ Mpc/h$'
        ylabel = r'$r_p\ /\ Mpc/h$'

        corr_object.plot_grid('',r_power,vmax=vmax,xlabel=xlabel,ylabel=ylabel)

    return

def bins_from_boundaries(boundaries):

    bins = list(zip(boundaries[:-1],boundaries[1:]))

    return bins

def plot_wedges(fig,ax,plot_info):

    #Unpack the plot_info dictionary.
    corr_obj = plot_info['corr_object']
    artists = []
    labels = []

    for i,mubin in enumerate(plot_info['mu_bins']):
        #Plot the data.
        colour = plot_info['mu_bin_colours'][i]
        try:
            abs_mu = plot_info['abs_mu']
        except:
            abs_mu = False
        if abs_mu:
            plot_label = r'${}<|\mu|<{}$'.format(mubin[0],mubin[1])
        else:
            plot_label = r'${}<\mu<{}$'.format(mubin[0],mubin[1])
        ar,lab = corr_obj.plot_wedge(ax,mubin,plot_label,colour,**plot_info['plot_data'],abs_mu=abs_mu)
        artists += [ar]
        labels += [lab]

        #Add a model or fit.
        if plot_info['plot_picca_fit']:
            rmin_fit = plot_info['picca_fit_data']['rmin']
            rmax_fit = plot_info['picca_fit_data']['rmax']
            try:
                fit_plot_data = plot_info['fit_plot_data']
            except:
                fit_plot_data = plot_info['plot_data']
            corr_obj.plot_wedge_fit(ax,mubin,'',colour,**fit_plot_data,rmin_fit=rmin_fit,rmax_fit=rmax_fit,abs_mu=abs_mu)
        if plot_info['plot_manual_fit']:
            b1,b2,beta1,beta2 = plot_info['manual_fit_data'].values()
            try:
                fit_plot_data = plot_info['fit_plot_data']
            except:
                fit_plot_data = plot_info['plot_data']
            corr_obj.plot_wedge_manual_model(ax,b1,b2,beta1,beta2,mubin,'',colour,fit_plot_data)

    #Add axis labels.
    if plot_info['format']['xlabel']:
        ax.set_xlabel(r'$r\ [\mathrm{{Mpc}}/h]$')
    if plot_info['format']['ylabel']:
        if plot_info['plot_data']['r_power'] > 1:
            ax.set_ylabel(r'$r^{:d}\ \xi (r)\ [(\mathrm{{Mpc}}\ h^{{-1}})^{{{:d}}}]$'.format(int(plot_info['plot_data']['r_power']),int(plot_info['plot_data']['r_power'])))
        elif plot_info['plot_data']['r_power'] == 1:
            ax.set_ylabel(r'$r\ \xi (r)\ [\mathrm{{Mpc}}\ h^{{-1}}]$')
        elif plot_info['plot_data']['r_power'] == 0:
            ax.set_ylabel(r'$\xi (r)$')

    #Add a legend if desired.
    if plot_info['format']['legend']:
        try:
            leg_loc = plot_info['format']['leg_loc']
        except:
            leg_loc = 0
        if leg_loc == 'shared':
            fig.legend(artists,labels,loc='lower center',borderaxespad=0,bbox_to_anchor=(0.5,0.05),ncol=len(plot_info['mu_bins']))
        else:
            ax.legend(loc=leg_loc)

    #Add a title if desired.
    if 'title' in list(plot_info['format'].keys()):
        if plot_info['format']['title'] is not None:
            ax.set_title(plot_info['format']['title'])

    return

def plot_rp_bins_vs_rt(ax,plot_info):

    #Unpack the plot_info dictionary.
    corr_obj = plot_info['corr_object']

    for i,rpbin in enumerate(plot_info['rp_bins']):
        #Plot the data.
        plot_label = r'${}<r_\parallel<{}$'.format(rpbin[0],rpbin[1])
        colour = plot_info['rp_bin_colours'][i]
        corr_obj.plot_rp_bin_vs_rt(ax,rpbin,plot_label,colour,**plot_info['plot_data'])

        #Add a model or fit.
        if plot_info['plot_picca_fit']:
            corr_obj.plot_rp_bin_vs_rt_fit(ax,rpbin,'',colour,**plot_info['plot_data'])
        if plot_info['plot_manual_fit']:
            b1,b2,beta1,beta2 = plot_info['manual_fit_data'].values()
            corr_obj.plot_rp_bin_vs_rt_manual_model(ax,b1,b2,beta1,beta2,rpbin,'',colour,**plot_info['plot_data'])

    #Add axis labels.
    if plot_info['format']['xlabel']:
        ax.set_xlabel(r'$r_\perp\ [\mathrm{{Mpc}}/h]$')
    if plot_info['format']['ylabel']:
        if plot_info['plot_data']['r_power'] > 1:
            ax.set_ylabel(r'$r^{:d}\ \xi (r)\ [(\mathrm{{Mpc}}\ h^{{-1}})^{{{:d}}}]$'.format(int(plot_info['plot_data']['r_power']),int(plot_info['plot_data']['r_power'])))
        elif plot_info['plot_data']['r_power'] == 1:
            ax.set_ylabel(r'$r\ \xi (r)\ [\mathrm{{Mpc}}\ h^{{-1}}]$')
        elif plot_info['plot_data']['r_power'] == 0:
            ax.set_ylabel(r'$\xi (r)$')

    #Add a legend if desired.
    if plot_info['format']['legend']:
        ax.legend()

    #Add a title if desired.
    if plot_info['format']['title'] is not None:
        ax.set_title(plot_info['format']['title'])

    return

def plot_rt_bins_vs_rp(ax,plot_info):

    #Unpack the plot_info dictionary.
    corr_obj = plot_info['corr_object']

    for i,rtbin in enumerate(plot_info['rt_bins']):
        #Plot the data.
        plot_label = r'${}<r_\parallel<{}$'.format(rtbin[0],rtbin[1])
        colour = plot_info['rt_bin_colours'][i]
        corr_obj.plot_rt_bin_vs_rp(ax,rtbin,plot_label,colour,**plot_info['plot_data'])

        #Add a model or fit.
        if plot_info['plot_picca_fit']:
            corr_obj.plot_rt_bin_vs_rp_fit(ax,rtbin,'',colour,**plot_info['plot_data'])
        if plot_info['plot_manual_fit']:
            b1,b2,beta1,beta2 = plot_info['manual_fit_data'].values()
            corr_obj.plot_rt_bin_vs_rp_manual_model(ax,b1,b2,beta1,beta2,rtbin,'',colour,**plot_info['plot_data'])

    #Add axis labels.
    if plot_info['format']['xlabel']:
        ax.set_xlabel(r'$r_\parallel\ [\mathrm{{Mpc}}/h]$')
    if plot_info['format']['ylabel']:
        if plot_info['plot_data']['r_power'] > 1:
            ax.set_ylabel(r'$r^{:d}\ \xi (r)\ [(\mathrm{{Mpc}}\ h^{{-1}})^{{{:d}}}]$'.format(int(plot_info['plot_data']['r_power']),int(plot_info['plot_data']['r_power'])))
        elif plot_info['plot_data']['r_power'] == 1:
            ax.set_ylabel(r'$r\ \xi (r)\ [\mathrm{{Mpc}}\ h^{{-1}}]$')
        elif plot_info['plot_data']['r_power'] == 0:
            ax.set_ylabel(r'$\xi (r)$')

    #Add a legend if desired.
    if plot_info['format']['legend']:
        ax.legend()

    #Add a title if desired.
    if plot_info['format']['title'] is not None:
        ax.set_title(plot_info['format']['title'])

    return
