import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from picca import wedgize
import sys
from numpy import linalg
import mcfit
from . import correlation_model


# TODO: write this
class picca_correlation:
    def __init__(parameters,correlation_data):

        #picca parameters
        self.correl_type = correl_type
        self.N_side = N_side
        self.N_pixels = N_pixels
        self.quantity_1 = quantity_1
        self.quantity_2 = quantity_2
        self.rpmax = rpmax
        self.rpmin = rpmin
        self.rtmax = rtmax
        self.np = np
        self.nt = nt
        self.rmax = rmax
        self.nr = nr
        self.zmin = zmin
        self.zmax = zmax

        #input skewer parameters
        self.sr = sr
        self.bm = bm

        #correlation data
        self.rp = rp
        self.rt = rt
        self.z = z
        self.xi_grid = xi_grid
        self.cov_grid = cov_grid

        #cosmology data
        self.Z_cosmology = Z_cosmology
        self.D_cosmology = D_cosmology

        #xcf_exp_1000_GAUSS_sr2.0_quantityGq_RN_nside16_rpmin-60.0_rpmax60.0_rtmax60.0_np40_nt20_zmin2.0_zmax2.2_bm1_biasG18_picos.fits.gz
        return

    @classmethod
    def make_correlation_object(cls,location,filename):

        #get parameters
        #replace with get_parameters_from_param_file when written
        parameters = get_parameters_from_filename(location+filename)

        #get correlation data
        h = fits.open(file_path)
        correlation_data = h[1].data

        #get cosmology data


        return cls(parameters,correlation_data)

    def plot_xir2(self,mubin,plot_label):

        mumin = mubin[0]
        mumax = mubin[1]

        b = wedgize.wedge(mumin=mumin,mumax=mumax,rtmax=self.rtmax,nrt=self.nt,
            rpmin=self.rpmin,rpmax=self.rpmax,nrp=self.np,nr=self.nr,
            rmax=self.rmax)

        #REMOVE NON_DIAGONALS FROM COV MATRIX
        for i in range(self.cov_grid.shape[0]):
            for j in range(self.cov_grid.shape[1]):
                if i!=j:
                    self.cov_grid[i,j]=0

        r, xi_wed, cov_wed = b.wedge(self.xi_grid,self.cov_grid)

        Nr = len(r)
        err_wed = np.zeros(Nr)
        for i in range(Nr):
            err_wed[i] = np.sqrt(cov_wed[i][i])

        cut = err_wed>0

        plt.errorbar(r[cut],xi_wed[cut],yerr=err_wed[cut],fmt='o',label=plot_label)

        return



def get_parameters_from_param_file(location,filename):

    # TODO: write this
    # TODO: in the picca function, make a param file too (use colore method for inspiration)

    return



def bins_from_boundaries(boundaries):

    bins = []
    for i in range(len(boundaries)-1):
        bins += [(boundaries[i],boundaries[i+1])]

    return bins

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
                r_CAMB,xi_CAMB_plot_label_CAMB = get_CAMB_xi(CAMB_location+CAMB_filename_format.format(sr),scale_CAMB[i],sr)
                plt.plot(r_CAMB,xi_CAMB*(r_CAMB**2),label=plot_label_CAMB)

        #save figure
        plt.legend(loc=3)
        plt.savefig('xi_wedge_{}_{}.pdf'.format(mumin,mumax))

def plot_per_file(mubins,filenames,add_CAMB,plot_label_parameters,CAMB_sr=None,scale_CAMB=None,CAMB_location=None,CAMB_filename_format=None):

    for filename in filenames:

        plt.figure()
        plt.axhline(y=0,color='gray',ls=':')
        plt.xlabel('r [Mpc/h]')
        plt.ylabel('xi(r)')
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
                add_CAMB_xi(CAMB_location+CAMB_filename.format(sr),scale_CAMB[i],sr)

        #save figure
        plt.legend(loc=3)
        plt.savefig('xi_wedges_{}.pdf'.format(plot_title))
