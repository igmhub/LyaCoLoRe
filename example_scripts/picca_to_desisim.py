import numpy as np
import fitsio
import matplotlib.pyplot as plt
import lya_mock_functions as mock

# identify output file we want to plot
h = fitsio.FITS('../example_data/delta_picca/z1.85/z1.85_N1000_node_015_nside_4_pix_10.fits')

# get information about quasars (TYPE,RA,DEC,Z_COSMO,DZ_RSD)
catalog = h[3].read()
z_qso = catalog['Z']
Nq = len(z_qso)
print('# quasars =',Nq)
print(np.min(z_qso),'< z_qso <',np.max(z_qso))

# get arraw with redshift in each cell of grid
loglam = h[2].read()
z = loglam/1215.67-1.0
Nz=len(z)

# Get deltas (fluctuation around mean density) 
delta = h[0].read()

filename='test_desisim_lya_N'+str(Nq)+'.fits'
fits = fitsio.FITS(filename,'rw',clobber=True)

# write skewers for desisim
for i in range(Nq):
    # Convert density to flux
    tau = mock.get_tau(z,1+delta[:,i])
    flux = np.exp(-tau)
    # only add absorption in the forest 
    no_forest = (z > z_qso[i])
    flux[no_forest]=1.0
    data = {}
    data['LAMBDA']=1215.67*(1+z)
    data['FLUX']=flux
    head = {}
    head['ZQSO']=z_qso[i]
    head['RA']=catalog['RA'][i]
    head['DEC']=catalog['DEC'][i]
    head['MAG_G']=22
    fits.write(data,header=head)
fits.close()


