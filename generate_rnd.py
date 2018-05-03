import numpy as np
from desimodel.footprint import is_point_in_desi
from desimodel.io import load_tiles
import astropy.table
from scipy.interpolate import interp1d

def generate_rnd(factor=3, out_path= None):
    """
    Routine to generate a random catalog in 3D following
    certain N(z) distribution
    
    Args:
    ----
    rad: z-values of the data catalog
    factor: Size of the generated catalog (before masking)
    out_path: Name of output file where randoms will be saved (default: None)

    Returns:
    -------
    x,y,z: Cartesian coordinates of the random points
    """
    #Creating random that follows N(z)
    tiles = load_tiles()
    nz_file = astropy.table.Table.read('example_data/Nz_qso_2_highZ.txt',format='ascii')
    zvec = np.linspace(np.min(nz_file['col1']),np.max(nz_file['col1']),5000)
    spl_z = interp1d(nz_file['col1'],nz_file['col2'])
    dndz = spl_z(zvec)
    ntot = 4*180**2/np.pi*np.sum(nz_file['col2'])*(nz_file['col1'][1]-nz_file['col1'][0])
    print 'Generating', factor*ntot, 'random points'
    z_rnd = np.random.choice(zvec[zvec>1.8],p=dndz[zvec>1.8]/np.sum(dndz[zvec>1.8]),size=int(factor*ntot))
    ra_rnd = 360.*np.random.random(size=len(z_rnd))
    cth_rnd = -1+2.*np.random.random(size=len(z_rnd))
    dec_rnd = np.arcsin(cth_rnd)*180/np.pi
    good = is_point_in_desi(tiles,ra_rnd,dec_rnd)
    ra_rnd = ra_rnd[good]
    dec_rnd = dec_rnd[good]
    z_rnd = z_rnd[good]
    if out_path is not None:
        tab_out = astropy.table.Table([ra_rnd,dec_rnd,z_rnd],names=('RA','DEC','Z'))
        tab_out.write(out_path,overwrite=True)
    return ra_rnd, dec_rnd, z_rnd
_, _, _ = generate_rnd(factor=1,out_path='test.fits.gz')
