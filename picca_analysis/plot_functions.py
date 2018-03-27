import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from picca import wedgize
import sys
from numpy import linalg
import mcfit

def get_parameters_from_filename(filename):

    default_cf_parameters = {'correlation': 'cf', 'nside': 8, 'sr': 2.0, 'rpmax': 200.0, 'rpmin': 0.0, 'rtmax': 200.0, 'np': 50, 'nt': 50, 'nr': 100, 'rmax': 200.0, 'zmin': 0.0, 'zmax': 4.0, 'quantity': 'GG'}
    default_xcf_parameters = {'correlation': 'xcf', 'nside': 8, 'sr': 2.0, 'rpmax': 200.0, 'rpmin': -200.0, 'rtmax': 200.0, 'np': 100, 'nt': 50, 'nr': 100, 'rmax': 200.0, 'zmin': 0.0, 'zmax': 4.0, 'quantity': 'Gq'}

    first_underscore = filename.find('_')
    if filename[first_underscore-3:first_underscore] == 'xcf':
        parameters = default_xcf_parameters
    elif filename[first_underscore-2:first_underscore] == 'cf':
        parameters = default_cf_parameters

    for parameter in parameters.keys():
        if parameter in filename:
            i_pvs = int(filename.find(parameter)+len(parameter))
            i_pve = int(i_pvs + filename[i_pvs:].find('_'))
            if i_pve < i_pvs:
                i_pve = int(i_pvs + filename[i_pvs:].find('.fits'))
            d_type = type(parameter)

            parameters[parameter] = float(filename[i_pvs:i_pve])

    return parameters

def get_plot_data(mumin,mumax,filename):

    #Read data from the filename to know how to wedgize etc
    #filename format: ${CORRELTYPE}_${NPIXELS}_${BRANCH}_${OPTION}_${RSDOPTION}_rpmin${RPMIN}_rpmax${RPMAX}_rtmax${RTMAX}_np${NP}_nt${NT}_zmin${ZQSOMIN}_zmax${ZQSOMAX}.fits.gz
    parameters = get_parameters_from_filename(filename)

    #file_parameters = {**file_parameters,**{filename:parameters}}

    h = fits.open(filename)
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

    r, xi_wed, cov_wed = b.wedge(xi_grid,cov_grid)

    Nr = len(r)
    err_wed = np.zeros(Nr)
    for i in range(Nr):
        err_wed[i] = np.sqrt(cov_wed[i][i])

    cut = err_wed>0

    return r, xi_wed, err_wed, cut, parameters

def add_CAMB_xi(location,filename,scale_CAMB=1,z=0.0,quantities='GG'):
    data = np.loadtxt(location+filename)

    #Get from input files
    bias_data = np.loadtxt()
    z_b = bias_data[0,:]
    b = bias_data[1,:]
    #Get from master file
    growth_data = np.loadtxt()
    z_D = growth_data[0,:]
    D = growth_data[1,:]

    b_at_zval = np.interp(z,z_b,b)
    D_at_zval = np.interp(z,z_D,D)
    D_at_z0 = np.interp(0,z_D,D)

    if quantities == 'GG':
        scale_CAMB = 1
    elif quantities == 'Gq':
        scale_CAMB = b_at_zval*D_at_zval/D_at_z0
    elif quantities == 'Dq':
        scale_CAMB = b_at_zval*(D_at_zval/D_at_z0)**2
    else:
        error('quantities not recognised')

    xi = data[:,1]
    r = data[:,0]
    xir2 = xi*r**2
    xir2 = xir2*scale_CAMB
    plt.plot(r,xir2,label='CAMB')

    return

def add_CAMB_Pk(location,filename,z,mu,RSDOPTION='NO_RSD',CAMB_input='xi'):
    data = np.loadtxt(location+filename)
    k = data[:,0]
    Pk = data[:,1]
    return

def plot_per_bin(mubins,filenames,CAMB_sr):

    for mubin in mubins:

        plt.figure()
        mumin = mubin[0]
        mumax = mubin[1]

        for filename in filenames:

            r,xi_wed,err_wed,cut,parameters = get_plot_data(mumin,mumax,filename)

            plot_label = ''
            for parameter in plot_label_parameters:
                plot_label += '{}: {}, '.format(parameter,parameters[parameter])
            plot_label = plot_label[:-2]

            plt.errorbar(r[cut],xi_wed[cut]*(r[cut]**2),yerr=err_wed[cut]*(r[cut]**2),fmt='o',label=plot_label)

            plt.axhline(y=0,color='gray',ls=':')
            plt.xlabel('r [Mpc/h]')
            plt.ylabel('r^2 xi(r)')
            plt.grid(True, which='both')
            plt.title('correlation function, {} < mu < {}'.format(mumin,mumax))

        if add_CAMB == True:
            for i,sr in enumerate(CAMB_sr):
                print(scale_CAMB[i]);add_CAMB_to_plot(CAMB_location,CAMB_filename.format(sr),scale_CAMB[i])

        #save figure
        plt.legend(loc=3)
        plt.savefig('xi_wedge_{}_{}.pdf'.format(mumin,mumax))

        #plt.ylim(2,20)
        #plt.xlim(70,120)
        #plt.savefig('xi_wedge_{}_{}_BAO_zoom.pdf'.format(mumin,mumax))

def plot_per_file(mubins,filenames,CAMB_sr):

    for filename in filenames:

        plt.figure()

        for mubin in mubins:

            mumin = mubin[0]
            mumax = mubin[1]
            r,xi_wed,err_wed,cut,parameters = get_plot_data(mumin,mumax,filename)

            plt.errorbar(r[cut],xi_wed[cut]*(r[cut]**2),yerr=err_wed[cut]*(r[cut]**2),fmt='o',label='{} < mu < {}'.format(mumin,mumax))

            plt.axhline(y=0,color='gray',ls=':')
            plt.xlabel('r [Mpc/h]')
            plt.ylabel('r^2 xi(r)')
            plt.grid(True, which='both')
            plt.title(filename[11:-8])

        if add_CAMB == True:
            for i,sr in enumerate(CAMB_sr):
                add_CAMB_to_plot(CAMB_location,CAMB_filename.format(sr),scale_CAMB[i])

        #save figure
        plt.legend(loc=3)
        plt.savefig('plot_{}.pdf'.format(filename[11:-8]))
