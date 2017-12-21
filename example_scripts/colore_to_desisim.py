import numpy as np
import fitsio
import matplotlib.pyplot as plt

# identify output file we want to plot
h = fitsio.FITS('../example_data/raw_colore/N1000_out_srcs_s0_15.fits')

# get information about quasars (TYPE,RA,DEC,Z_COSMO,DZ_RSD)
catalog = h[1].read()
z_qso = catalog['Z_COSMO']
Nq = len(z_qso)
print('# quasars =',Nq)
print(np.min(z_qso),'< z_qso <',np.max(z_qso))

# get arraw with redshift in each cell of grid
z = h[4].read()['Z']
Nz=len(z)

# Get deltas (fluctuation around mean density) and line of sight velocity 
delta = h[2].read()
velocity = h[3].read()
# Convert density to flux
import lya_mock_functions as mock
tau = mock.get_tau(z,1+delta)
flux = np.exp(-tau)

filename='desisim_lya_N'+str(Nq)+'.fits'
fits = fitsio.FITS(filename,'rw',clobber=True)

# write skewers for desisim
for i in range(Nq):
    # only add absorption in the forest 
    no_forest = (z > z_qso[i])
    flux[i][no_forest]=1.0
    data = {}
    data['LAMBDA']=1215.67*(1+z)
    data['FLUX']=flux[i]
    head = {}
    head['ZQSO']=z_qso[i]
    head['RA']=catalog['RA'][i]
    head['DEC']=catalog['DEC'][i]
    head['MAG_G']=22
    fits.write(data,header=head)
fits.close()
