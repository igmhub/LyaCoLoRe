import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import fitsio
import healpy as hp
import os

from lyacolore import utils

# Set up colours
mycolours = {'C0': '#F5793A', 'C1': '#A95AA1', 'C2': '#85C0F9', 'C3': '#0F2080'}

# Plot settings
fontsize = 18
plt.rc('font', size=fontsize)

def quant_hist(cats,cattypes,labels=None,quant='z',density=True):
    
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
        
    for i,cat in enumerate(cats):
        
        cattype = cattypes[i]
        assert cattype in ['master','master_randoms','master_dla','master_dla_randoms','zcat','drq']
        
        h = fitsio.FITS(cat)
        if cattype=='master':
            z = h[1][:]['Z_QSO_RSD']
        elif cattype=='master_randoms':
            z = h[1][:]['Z']
        elif cattype=='master_dla':
            z = h[1][:]['Z_DLA_RSD']
        elif cattype=='master_dla_randoms':
            z = h[1][:]['Z_DLA']
        elif cattype=='zcat':
            z = h[1][:]['Z']
            gm = 22.5 - 2.5 * np.log10(h[1][:]['FLUX_G'])
            rm = 22.5 - 2.5 * np.log10(h[1][:]['FLUX_R'])
            zm = 22.5 - 2.5 * np.log10(h[1][:]['FLUX_Z'])
        elif cattype=='drq':
            z = h[1][:]['Z']
        h.close()

        if labels is not None:
            label = labels[i]
        else:
            label = None
            
        if quant == 'z':
            q = z
            m, M = 1.5, 4.0
            qname = 'z'
        elif quant == 'gmag':
            assert cattype == 'zcat'
            q = gm
            m, M = 18., 23.
            qname = 'm_G'
        elif quant == 'rmag':
            assert cattype == 'zcat'
            q = rm
            m, M = 18., 23.
            qname = 'm_R'
        elif quant == 'zmag':
            assert cattype == 'zcat'
            q = zm
            m, M = 18., 23.
            qname = 'm_Z'
        qname = r'$\langle {} \rangle$'.format(qname)
                    
        bins = np.linspace(m,M,501)
        ax.hist(q,bins=bins,histtype='step',density=density,label=label,color=mycolours['C{}'.format(i)])
        
    ax.set_xlabel(qname)
    ax.legend()
    
    return
    

def mean_quant_sky(cat,cattype,nside=64,quant='z',vmin=2.0,vmax=2.5,zmin=None,title=None):
    
    ## Get key healpix quantities
    npix = hp.nside2npix(nside)
    pixarea = hp.nside2pixarea(nside,degrees=True)
    
    ## Check that the args are ok
    assert cattype in ['master','master_randoms','zcat','drq']
    assert quant in ['z','ndens','gmag','rmag','zmag']
    
    ## Open the catalogue and extract data
    h = fitsio.FITS(cat)
    if cattype=='master':
        ra = h[1][:]['RA']
        dec = h[1][:]['DEC']
        z = h[1][:]['Z_QSO_RSD']
    elif cattype=='master_randoms':
        ra = h[1][:]['RA']
        dec = h[1][:]['DEC']
        z = h[1][:]['Z']
    elif cattype=='zcat':
        ra = h[1][:]['RA']
        dec = h[1][:]['DEC']
        z = h[1][:]['Z']
        gm = 22.5 - 2.5 * np.log10(h[1][:]['FLUX_G'])
        rm = 22.5 - 2.5 * np.log10(h[1][:]['FLUX_R'])
        zm = 22.5 - 2.5 * np.log10(h[1][:]['FLUX_Z'])
    elif cattype=='drq':
        ra = h[1][:]['RA']
        dec = h[1][:]['DEC']
        z = h[1][:]['Z']
    h.close()
    
    ## Apply the min redshift if needed
    if zmin is not None:
        w = z>=zmin
        ra = ra[w]
        dec = dec[w]
        z = z[w]
        if cattype=='zcat':
            gm = gm[w]
            rm = rm[w]
            zm = zm[w]
    
    ## Compute the healpix pixels
    theta = (90-dec)*np.pi/180.
    phi = ra*np.pi/180.
    hpix = hp.ang2pix(nside,theta,phi,nest=False)

    ## Construct the vector for our output quantity
    q = np.zeros(12*nside*nside)
    counts = np.bincount(hpix,minlength=npix)
    w = counts>0

    ## Compute the output quantity
    if quant == 'z':
        zsums = np.bincount(hpix,weights=z,minlength=npix)
        q[w] = zsums[w]/counts[w]
        qname = r'$\langle z \rangle$'
    elif quant == 'ndens':
        q[w] = counts[w]/pixarea
        qname = r'$\langle n \rangle~[\mathrm{deg}^{-2}]$'
    elif quant == 'gmag':
        gmsums = np.bincount(hpix,weights=gm,minlength=npix)
        q[w] = gmsums[w]/counts[w]
        qname = r'$\langle m_G \rangle$'
    elif quant == 'rmag':
        rmsums = np.bincount(hpix,weights=rm,minlength=npix)
        q[w] = rmsums[w]/counts[w]
        qname = r'$\langle m_R \rangle$'
    elif quant == 'zmag':
        zmsums = np.bincount(hpix,weights=zm,minlength=npix)
        q[w] = zmsums[w]/counts[w]
        qname = r'$\langle m_Z \rangle$'
            
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111,projection='mollweide')
    
    phi_grid = np.linspace(0,2*np.pi,1000)
    theta_grid = np.linspace(0,np.pi,1000)
    lon_grid = phi_grid - np.pi
    lat_grid = theta_grid - np.pi/2.
    
    PHI, THETA = np.meshgrid(phi_grid, theta_grid)
    grid_pix = hp.ang2pix(nside, THETA, PHI)
    
    grid_map = np.ma.array(q[grid_pix],mask=(~w)[grid_pix])
    cmap =  matplotlib.cm.get_cmap('plasma')
    image = plt.pcolormesh(lon_grid[::-1], lat_grid[::-1], grid_map, vmin=vmin, vmax=vmax, rasterized=True, cmap=cmap)
        
    ax.set_longitude_grid(60)
    
    cb = fig.colorbar(image, orientation='horizontal', shrink=.6, pad=0.05, ticks=[vmin, vmax])
    cb.ax.xaxis.set_label_text(qname)
    cb.ax.xaxis.labelpad = -8
    
    ax.grid()
    
    if title is not None:
        ax.set_title(title)
    
    plt.tight_layout()
    plt.show()
    
    return