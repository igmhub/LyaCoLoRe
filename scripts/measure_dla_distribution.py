import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

from lyacolore import DLA

#Import pyigm modules
from pyigm.fN.fnmodel import FNModel
fN_default = FNModel.default_model()
fN_default.zmnx = (0.5,5)
fN_cosmo = fN_default.cosmo
use_pyigm = True

include_RSDs = True

#Open a DLA master file and extract the relevant data.
h = fits.open(basedir+'/master_DLA.fits')
if include_RSDs:
    dla_z = h['DLACAT'].data['Z_DLA_RSD']
else:
    dla_z = h['DLACAT'].data['Z_DLA_NO_RSD']
dla_NHI = h['DLACAT'].data['N_HI_DLA']
h.close()

#Open the master file to get QSO information (needed for determining absorption length)
h = fits.open(basedir+'/master.fits')
z_qso = h['CATALOG'].data['Z_QSO_NO_RSD']
z_cell = h['COSMO_COl'].data['Z']
N_qso = z_qso.shape[0]
h.close()

#Function to measure dn/dz
def measure_dndz(dla_z,dla_NHI,log_NHI_bin,z_bins):

    #Filter the DLAs according to the bin in log NHI.
    log_NHI_min = log_NHI_bin[0]
    log_NHI_max = log_NHI_bin[1]
    relevant_dlas = (dla_NHI>log_NHI_min) * (dla_NHI<log_NHI_max)
    dla_z = dla_z[relevant_dlas]

    n = np.zeros(len(z_bins))
    z = np.zeros(len(z_bins))
    z_widths = np.zeros(len(z_bins))

    for k,z_bin in enumerate(z_bins):
        z_min = z_bin[0]
        z_max = z_bin[1]

        #Use the QSO data to determine the effective redshift interval over all QSOs.
        z_low = np.ones(z_qso.shape)*z_min
        z_high = np.minimum(np.maximum(z_qso,z_max),z_max)
        dz_total = np.sum(z_high-z_low)
        z_widths[k] = dz_total

        #Use the QSO data to get the average redshift of the bin.
        N_qso_in_z_bin = np.sum((z_qso<=z_max)*(z_qso>z_min))
        z[k] = z_min + dz_total/(2*N_qso_in_z_bin)

        #Count the DLAs in the bin.
        dlas_in_z_bin = (dla_z<=z_max)*(dla_z>z_min)
        n[k] = np.sum(dlas_in_z_bin)

    #Calculate dndz
    dndz = n/(z_widths)

    return z,dndz

#Function to measure f(NHI)
def measure_fNHI(dla_z,dla_NHI,z_bin,log_NHI_bin_edges):

    #Make the NHI bins.
    log_NHI_bins = []
    for i in range(log_NHI_bin_edges.shape[0]-1):
        log_NHI_bins += [(log_NHI_bin_edges[i],log_NHI_bin_edges[i+1])]
    log_NHI_bins = np.array(log_NHI_bins)
    NHI_bin_edges = 10**(log_NHI_bin_edges)
    NHI_bin_widths = NHI_bin_edges[1:] - NHI_bin_edges[:-1]

    #Filter the DLAs by z bin
    z_min = z_bin[0]
    z_max = z_bin[1]
    dlas_in_z_bin = (dla_z<=z_max)*(dla_z>=z_min)

    #Use the QSO data to determine the effective redshift interval over all QSOs.
    z_low = np.ones(z_qso.shape)*z_min
    z_high = np.minimum(np.maximum(z_qso,z_max),z_max)
    dz_total = np.sum(z_high-z_low)

    #Use the QSO data to get the average redshift of the bin.
    N_qso_in_z_bin = np.sum((z_qso<=z_max)*(z_qso>z_min))
    z_val = z_min + dz_total/(2*N_qso_in_z_bin)

    #Calculate dX
    dX = fN_cosmo.abs_distance_integrand(z_val)  *  dz_total

    #In each of the NHI bins, count the DLAs
    n_values = []
    log_NHI_values = []

    for i,lNb in enumerate(log_NHI_bins):
        log_NHI_min = lNb[0]
        log_NHI_max = lNb[1]
        counters = (dla_NHI[dlas_in_z_bin]>log_NHI_min)*(dla_NHI[dlas_in_z_bin]<log_NHI_max)
        n_value = np.sum(counters)
        n_values += [n_value]
        if n_value>0:
            log_NHI_eff = np.log10(np.average(10**dla_NHI[dlas_in_z_bin][counters]))
        else:
            log_NHI_eff = (log_NHI_max+log_NHI_min)/2.
        log_NHI_values += [log_NHI_eff]

    n_values = np.array(n_values)
    log_NHI_values = np.array(log_NHI_values)

    #Calculate f by dividing by the bin width and the effective absorption length (dX in each of Nq skewers)
    f_values = n_values/(Nq*NHI_bin_widths*dX)

    return z_val,log_NHI_values,f_values


#Measure dndz and check that it matches the input.
z_edges = np.linspace(1.5,4.0,22)
z_bins = []
for i in range(z_edges.shape[0]-1):
    z_bins += [(z_edges[i],z_edges[i+1])]

plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')

N_intervals = [(17.2,22.0),(17.2,20.3),(20.3,22.0)]
for k,N_int in enumerate(N_intervals):

    z_values,dndz_measured = measure_dndz(dla_z,dla_NHI,N_int,z_bins)
    plt.scatter(z_values,dndz_measured,marker='x',label=r'${}<\log N_{{HI}}<{}$'.format(Nmin_bin,Nmax_bin))

    dndz_model = DLA.dndz(z_values,Nmin_bin,Nmax_bin)
    plt.plot(z_values,dndz_model,c='C{}'.format(k),linestyle='--')

plt.xlim(2.2,3.6)
plt.ylim(0.0,3.5)
plt.xlabel(r'$z$')
plt.ylabel(r'$\mathrm{{d}}n / \mathrm{{d}}z$')
plt.legend()
plt.grid()
plt.savefig('dndz.pdf')
plt.show()


#Measure f(NHI)
z_bins = [(2.1,2.3),(2.5,2.7),(2.9,3.1),(3.3,3.5)]
log_NHI_bin_edges = np.linspace(17.2,22.5,30)

plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')

for k,z_bin in enumerate(z_bins):

    #Measure the f values and scatter plot.
    z_val,log_NHI_values,f_values = measure_fNHI(dla_z,dla_NHI,z_bin,log_NHI_bin_edges,z_cell)
    non_empty = (f_values > 0)
    plt.scatter(log_NHI_values[non_empty],np.log10(f_values[non_empty]),label=r'$z={:1.2f}$'.format(z_val),c='C{}'.format(k),marker='x')

    #Compute the input from pyigm and plot a line,.
    log10_f_model = fN_default.evaluate(log_NHI_values,z_val).reshape(log_NHI_values.shape[0])
    plt.plot(log_NHI_values,log10_f_model,c='C{}'.format(k))

    #plt.plot(log_NHI_values,(f_values-10**log10_f_model)/(10**log10_f_model),label=r'$z={:1.2f}$'.format(z_val),c='C{}'.format(k))

#Control the axes etc of the plot.
plt.xlim(17.0,23.0)
plt.ylim(-28.0,-17.0)
#plt.xlim(18.5,20.0)
#plt.ylim(-21.0,-20.0)
#plt.ylim(-0.1,0.1)
plt.ylabel(r'$\log(f(N_{{HI}}))$')
plt.xlabel(r'$\log(N_{{HI}})$')
plt.legend()
plt.grid()
#plt.savefig('f_NHI.pdf')
plt.show()
