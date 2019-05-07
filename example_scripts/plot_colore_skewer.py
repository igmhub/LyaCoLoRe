import numpy as np
import fitsio
import matplotlib.pyplot as plt
import lya_mock_functions as mock

# identify output file we want to plot
h = fitsio.FITS('../example_data/raw_colore_1000/out_srcs_s1_0.fits')

# get information about quasars (TYPE,RA,DEC,Z_COSMO,DZ_RSD)
catalog = h[1].read()
z_qso = catalog['Z_COSMO']
Nq = len(z_qso)
print('# quasars =',Nq)
print(np.min(z_qso),'< z_qso <',np.max(z_qso))

# get arraw with redshift in each cell of grid
z = h[4].read()['Z']
Nz=len(z)

# get deltas (fluctuation around mean density) and line of sight velocity
delta = h[2].read()
velocity = h[3].read()

# get rid of low-z quasars that would make uninteresting plots
highz = np.where(z_qso > 2.5)
z_qso = z_qso[highz]
catalog = catalog[highz]
delta = delta[highz]
velocity = velocity[highz]

# plot density vs redshift
plt.plot(z,delta[0])
plt.xlabel('z')
plt.ylabel(r'$\delta$')
plt.title('Density fluctuations in CoLoRe skewer')
plt.show()

# plot line of sight velocity vs redshift
plt.plot(z,velocity[0])
plt.xlabel('z')
plt.ylabel('velocity')
plt.title('Line of sight velocity in CoLoRe skewer')
plt.show()

