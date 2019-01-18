import numpy as np
import astropy.table
import fitsio
import os
def generate_rnd(factor=3, out_path=None):
    """
    Routine to generate a random catalog in 3D following
    certain N(z) distribution
    
    Args:
    ----
    factor: Size of the generated catalog (before masking)
    out_path: Output path
    """
    #Creating random that follows N(z)
    data = fitsio.read('/global/projecta/projectdirs/desi/mocks/lya_forest/london/v5.0.0/master_DLA.fits')
    zvec = data['Z_DLA_RSD']
    ntot = int(len(data)*factor)
    z_rnd = np.random.choice(zvec,size=ntot)+0.0025*np.random.normal(size=ntot)
    qso_rnd = np.random.choice(np.arange(len(data)), size=ntot)
    ra_rnd = data['RA'][qso_rnd]
    dec_rnd = data['DEC'][qso_rnd]
    if out_path is not None:
        tab_out = astropy.table.Table([ra_rnd,dec_rnd,z_rnd],names=('RA','DEC','Z'))
        tab_out.write(out_path,overwrite=True)
    return None
# Execute
generate_rnd(factor=10,out_path='/global/projecta/projectdirs/desi/mocks/lya_forest/london/v5.0.0/master_DLA_randoms.fits.gz')
