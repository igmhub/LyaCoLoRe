import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import lya_mock_functions as mock
import lya_stats_functions as stats

#Open data file
hdulist = fits.open('/Users/jfarr/Projects/repixelise/test_input/cut_10000_out_srcs_s0_15.fits')

#Extract redshift from data file
z = hdulist[4].data['Z']
z = np.asarray(z)

#Get number of quasars, and redshift array
z_qso = hdulist[1].data['Z_COSMO']
N_qso = len(z_qso)
print('There are %d quasars in the sample.' % N_qso)

#Plot the locations of the quasars according to the HEALPix pixel grid.
RA = hdulist[1].data['RA']
DEC = hdulist[1].data['DEC']
phi = RA*np.pi/180
theta = np.pi/2 - DEC*np.pi/180
plt.figure()
plt.scatter(phi/np.pi,np.cos(theta))
plt.xlim(0.0,2.0)
plt.ylim(-1.0,1.0)

#Extract the delta skewers from the file, and make a mask for them.
delta_skewers = hdulist[2].data
mask = stats.make_mask(z,z_qso)

#NORMALISE THE DELTA FIELD SO THAT ITS MEAN IS 0.
#THIS IS AN ONGOING ISSUE IN COLORE, AND IS LISTED TO BE FIXED SO THAT THE DELTA MEAN IS 0 AUTOMATICALLY.
delta_skewers = stats.normalise_delta(mask,delta_skewers)

#Get the length of each skewer.
N_pix_skewer = delta_skewers.shape[1]
print('There are %d pixels in each skewer.' % N_pix_skewer)

#Show the structure of the data
print(hdulist[2].header.keys)

#Get delta for the highest redshift quasar in the sample.
id = np.argmax(hdulist[1].data['Z_COSMO'])
delta = delta_skewers[id,:]

#Show delta vs r
print('mean delta =', np.average(delta,weights=mask[id]))
print('var delta =', np.var(delta))
plt.figure()
plt.plot(hdulist[4].data['R'],delta)
plt.xlabel('$r\\,\\,[{\\rm Mpc}/h]$')
plt.ylabel('$\\delta$')

#Show delta vs z
plt.figure()
plt.plot(z,delta)
plt.xlabel('z')
plt.ylabel('$\\delta$')
print(" ")

#Show density vs z
density = (delta + 1)
print('mean density =', np.mean(density))
print('var density =', np.var(density))
plt.figure()
plt.semilogy(z,density)
plt.xlabel('z')
plt.ylabel('density')
print(" ")

#Convert from a lognormal density field to optical depth
tau = mock.get_tau(z,density)
print('mean tau =', np.mean(tau))
print('var tau =', np.var(tau))
plt.figure()
plt.semilogy(z,tau)
plt.xlabel('z')
plt.ylabel('optical depth')
print(" ")

#Convert from optical depth to transmitted flux fraction
flux = np.exp(-tau)
print('mean flux =', np.mean(flux))
print('var flux =', np.var(flux))
plt.figure()
plt.plot(z,flux)
plt.xlabel('z')
plt.ylabel('transmitted flux fraction')
plt.xlim(2,1.01*max(z))
print(" ")

#Show a test skewer spectrum
wavelength = 1215.67*(1+z)
plt.figure()
plt.plot(wavelength,flux)
plt.xlabel('wavelength')
plt.ylabel('transmitted flux fraction')
plt.xlim(1215.67*(1+2),1.01*max(wavelength))
print(" ")

#Calculate density statistics
binned_z, mean_binned_density, binned_density_var, binned_delta_var = stats.density_stats(z,z_qso,(delta_skewers)+1)

#Show the calculated statistics against z
plt.figure()
plot_mean_binned_density = plt.plot(binned_z, mean_binned_density)
plt.xlabel('z')
plt.ylabel('Mean density')
plt.figure()
plot_binned_density_var = plt.plot(binned_z, binned_density_var)
plt.xlabel('z')
plt.ylabel('Density variance')
plt.figure()
plot_binned_delta_var = plt.plot(binned_z, binned_delta_var)
plt.xlabel('z')
plt.ylabel('Delta variance')
print(" ")


#Show a sequence of histograms of delta
#Define redshift bin boundaries: zhb(0)<=z<zhb(1), zhb(1)<=z<zhb(2) etc.
z_hist_bins_boundaries = [0,2,3]
N_bins = len(z_hist_bins_boundaries)

lower_bound = [0]*N_bins
upper_bound = [0]*N_bins

for i in range(N_bins):
    #get boundaries in skewer length terms
    lower_bound[i] = max(np.argmax(z>z_hist_bins_boundaries[i])-1,0)
    if i+1 < N_bins:
        upper_bound[i] = np.argmax(z>z_hist_bins_boundaries[i+1])-1
    else:
        upper_bound[i] = len(z)-1

    #print histogram
    plt.figure()
    plt.hist(np.ravel(delta_skewers[:,lower_bound[i]:upper_bound[i]]),bins=1000,weights=np.ravel(mask[:,lower_bound[i]:upper_bound[i]]))
    plt.yscale('log',nonposy='clip')
    plt.xlabel('$\\delta$')
    plt.ylabel('frequency')

    if i+1 < N_bins:
        plt.title('{}<z<{}'.format(z_hist_bins_boundaries[i], z_hist_bins_boundaries[i+1]))
    else:
        plt.title('z>{}'.format(z_hist_bins_boundaries[i]))

#Make one plot with all of the
#plt.figure()
#plt.yscale('log',nonposy='clip')
#plt.xscale('log',nonposy='clip')
#for i in range(N_bins):
#    plt.hist(np.ravel(delta_skewers[:,lower_bound[i]:upper_bound[i]]),bins=100,weights=np.ravel(mask[:,lower_bound[i]:upper_bound[i]]),normed=True,histtype='step')
#plt.xlabel('$\\delta$')
#plt.ylabel('frequency')
#plt.xlim(-1,50)

plt.show()


#Old stuff below

#First HDU contains the source catalog
#print(hdulist[1].header.keys)
#plt.figure(); plt.hist(hdulist[1].data['Z_COSMO'],bins=100)
#print(" ")

#Second HDU contains the density skewers as a FITS image
#The skewers have the same ordering as the sources in the catalog
#(i.e. skewer delta_skewers[i,:] corresponds to source hdulist[1].data[i])

#Third HDU contains the velocity skewers. The units of the velocity are
#such that the skewers contain the redshift distortion associated with
#the peculiar velocity field
#print(hdulist[3].header.keys)
#plt.figure(); plt.plot(hdulist[4].data['R'],hdulist[3].data[id]);
#plt.xlabel('$r\\,\\,[{\\rm Mpc}/h]$',fontsize=18);
#plt.ylabel('$\\delta z_{\\rm RSD}$',fontsize=18)
#print(" ")

#Fourth HDU is a table containing background cosmological quantities at
#the distances where the skewers are sampled (see the use of
#hdulist[4].data['R'] in the previous examples
#print(hdulist[4].header.keys)
#plt.figure();
#plt.plot(hdulist[4].data['Z'],hdulist[4].data['R']*0.001,label='$r(z)\\,[{\\rm Gpc}/h]$')
#plt.plot(hdulist[4].data['Z'],hdulist[4].data['D'],label='$D_\\delta(z)$')
#plt.plot(hdulist[4].data['Z'],hdulist[4].data['V'],label='$D_v(z)$')
#plt.legend(loc='lower right')
#plt.xlabel('$z$',fontsize=18)
#print(" ")

#
