import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from picca import wedgize
import sys
from numpy import linalg

default_filenames = ['cf_exp.fits.gz']

#Set up the list of files to wedgize.
N_files = len(sys.argv) - 1
if N_files > 0:
    filenames = sys.argv[1:]
else:
    filenames = default_filenames

print('The cf_exp_out file(s) to wedgize are:')
for filename in filenames:
    print(filename)

#Set the default parameter values for wedgizing. These match picca's defaults
default_cf_parameters = {'nside': 8, 'sr': 2.0, 'rpmax': 200.0, 'rpmin': 0.0, 'rtmax': 200.0, 'np': 50, 'nt': 50, 'nr': 100, 'rmax': 200.0, 'zmin': 0.0, 'zmax': 4.0}
default_xcf_parameters = {'nside': 8, 'sr': 2.0, 'rpmax': 200.0, 'rpmin': -200.0, 'rtmax': 200.0, 'np': 100, 'nt': 50, 'nr': 100, 'rmax': 200.0, 'zmin': 0.0, 'zmax': 4.0}
file_parameters = {}

#Set the parameters to label plots by
plot_label_parameters = ['zmin','zmax']

#Either plot_per_file (i.e. each file gets its own plot, showing all mu bins)
#Or plot_per_bin (i.e. each mu bin has its own plot, but all files are grouped)
mode = 'plot_per_bin'

#Option to add in a scaled version of a CAMB power spectrum:
quantity = 'gaussian'
add_CAMB = False
CAMB_sr = ['10']
#CAMB_sr = ['10','10','10']
scale_CAMB = [1.0]
#scale_CAMB = [2.9488*0.4084/0.9998, 3.4069*0.3534/0.9998 , 3.8682*0.3110/0.9998]
#CAMB_location = '/Users/James/Projects/LyaCoLoRe/camb_scripts/'
CAMB_location = '/global/homes/j/jfarr/Projects/LyaCoLoRe/camb_scripts/'
CAMB_filename = 'camb_xi_{}.txt'

#Set up the bins of mu.
#ith bin is the range mubins[i]:mubins[i+1]
mubin_boundaries = [0.0,0.5,0.8,0.95,1.0]

mubins = []
for i in range(len(mubin_boundaries)-1):
    mubins += [(mubin_boundaries[i],mubin_boundaries[i+1])]
N_bins = len(mubins)

#Get the data from the wedge and filename to plot.
def get_plot_data(mumin,mumax,filename):

    #Read data from the filename to know how to wedgize etc
    #filename format: ${CORRELTYPE}_${NPIXELS}_${BRANCH}_${OPTION}_${RSDOPTION}_rpmin${RPMIN}_rpmax${RPMAX}_rtmax${RTMAX}_np${NP}_nt${NT}_zmin${ZQSOMIN}_zmax${ZQSOMAX}.fits.gz
    if filename[:filename.find('_')] == 'cf':
        parameters = default_cf_parameters
    elif filename[:filename.find('_')] == 'xcf':
        parameters = default_xcf_parameters

    for parameter in parameters.keys():
        if parameter in filename:
            i_pvs = int(filename.find(parameter)+len(parameter))
            i_pve = int(i_pvs + filename[i_pvs:].find('_'))
            if i_pve < i_pvs:
                i_pve = int(i_pvs + filename[i_pvs:].find('.fits'))
            d_type = type(parameter)

            parameters[parameter] = float(filename[i_pvs:i_pve])

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
    #w,v=linalg.eig(cov_grid)
    #print('min eig before non-diag removed',np.min(w))
    #REMOVE NON_DIAGONALS FROM COV MATRIX
    for i in range(cov_grid.shape[0]):
        for j in range(cov_grid.shape[1]):
            if i!=j:
                cov_grid[i,j]=0
    #w,v=linalg.eig(cov_grid)
    #print('min eig before non-diag removed',np.min(w))
    r,xi_wed,cov_wed = b.wedge(xi_grid,cov_grid)

    Nr = len(r)
    err_wed = np.zeros(Nr)
    for i in range(Nr):
        err_wed[i] = np.sqrt(cov_wed[i][i])
        #print(i,err_wed[i],cov_wed[i][i])

    cut = err_wed>0

    return r, xi_wed, err_wed, cut, parameters

def add_CAMB_to_plot(location,filename,scale_CAMB=1):
    data = np.loadtxt(location+filename)
    xi = data[:,1]
    r = data[:,0]
    xir2 = xi*r**2
    xir2 = xir2*scale_CAMB
    plt.plot(r,xir2,label='CAMB')
    return

#Make one plot per bin, with all files.
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

#Make one plot per file, with all bins.
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

if mode == 'plot_per_bin':
    plot_per_bin(mubins,filenames,CAMB_sr)
elif mode == 'plot_per_file':
    plot_per_file(mubins,filenames,CAMB_sr)
else:
    error('Mode not recognised.')

plt.show()
