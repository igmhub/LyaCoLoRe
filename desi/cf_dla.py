import numpy as np
import matplotlib.pyplot as plt
import treecorr as tc
import astropy.table
import healpy as hp
import glob
#from desimodel.footprint import is_point_in_desi
#from desimodel.io import load_tiles
import fitsio
# Cosmology
omega_matter = 0.140247/0.6800232**2
Omega_baryon = 0.022337/0.6800232**2
Omega_curvature = 0
H0 = 68.002320
sigma_8 = 0.811322
n_s = 0.963180
zmin_zmax = [(0,5.)]
from astropy.cosmology import FlatLambdaCDM
cosmo=FlatLambdaCDM(H0=H0,Om0=omega_matter)
path = '/global/projecta/projectdirs/desi/mocks/lya_forest/london/v4.0.0/DLA.fits.gz'
rnd_path = '/global/projecta/projectdirs/desi/mocks/lya_forest/london/v4.0.0/DLA_randoms.fits.gz'
source_table = fitsio.read(path)
redshift = source_table['ZDLA']
print 'Number of sources', len(source_table)
rad=cosmo.comoving_distance(redshift)*H0/100.
catx=rad*np.cos(source_table['DEC']/180*np.pi)*np.cos(source_table['RA']/180*np.pi)
caty=rad*np.cos(source_table['DEC']/180*np.pi)*np.sin(source_table['RA']/180*np.pi)
catz=rad*np.sin(source_table['DEC']/180*np.pi)
rnd = fitsio.read(rnd_path)
r_rnd = cosmo.comoving_distance(rnd['Z'])*H0/100.
ra_rnd = rnd['RA']
dec_rnd = rnd['DEC']
xrand=r_rnd*np.cos(ra_rnd*np.pi/180.)*np.cos(dec_rnd*np.pi/180.)
yrand=r_rnd*np.sin(ra_rnd*np.pi/180.)*np.cos(dec_rnd*np.pi/180.)
zrand=r_rnd*np.sin(dec_rnd*np.pi/180.)
icount=0
print(np.count_nonzero(np.isnan(xrand)), np.count_nonzero(np.isnan(yrand)), np.count_nonzero(np.isnan(zrand)))
for zmin,zmax in zmin_zmax:
    print "Doing",zmin,"-",zmax
    dimin,dimax=cosmo.comoving_distance([zmin,zmax])*H0/100.
    print dimin, dimax
    wh=(redshift>zmin) & (redshift<zmax)
    print 'Bin contains : ', np.count_nonzero(wh), ' galaxies'
    sigcat=tc.Catalog(x=catx[wh],y=caty[wh],z=catz[wh])
    whr = (rnd['Z']>zmin) & (rnd['Z']<zmax)
    print 'Bin contains : ', np.count_nonzero(whr), ' random objects'
    rancat=tc.Catalog(x=xrand[whr],y=yrand[whr],z=zrand[whr])
    dd=tc.NNCorrelation(min_sep=0.1,bin_size=0.075,max_sep=250.)
    dd.process(sigcat)
    dr=tc.NNCorrelation(min_sep=0.1,bin_size=0.075,max_sep=250.)
    dr.process(sigcat,rancat)
    rr=tc.NNCorrelation(min_sep=0.1,bin_size=0.075,max_sep=250.)
    rr.process(rancat,rancat)
    xi,xivar=dd.calculateXi(rr,dr)
    tab = astropy.table.Table([np.exp(dd.logr),xi,xivar],names=('r','xi','xivar'))
    filename = 'xi_dla_dla_v4.fits.gz'
    tab.write(filename,overwrite=True)
    icount=icount+1

