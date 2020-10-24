import numpy as np
import scipy as sp
import fitsio
import subprocess

basedir = "/project/projectdirs/desi/users/jfarr/LyaCoLoRe_paper/"
v_maj = 9
v_min = 0
v_realisations = range(10)
zbins = [(0.0,2.0),(2.0,2.2),(2.2,2.4),(2.4,2.6),(2.6,2.8),(2.8,3.0),(3.0,3.4),(3.4,10.0)]

# Function to coadd a set of cf of xcf measurements.
def coadd_correlations(fi,fout=None):
    # Define variables of the correct shape to store correlation information.
    h = fitsio.FITS(fi[0])
    head = h[1].read_header()
    nobj = 0
    rp = h[1]['RP'][:]*0
    rt = h[1]['RT'][:]*0
    nb = h[1]['NB'][:]*0
    z = h[1]['Z'][:]*0
    hid = h[2]['HEALPID'][:]
    wet = rp*0
    da = sp.zeros(h[2]['DA'][:].shape)
    we = sp.zeros(h[2]['WE'][:].shape)
    h.close()

    # For each data file:
    for f in fi:
        print("coadding file {}".format(f),end="\r")

        # Add information about the weights and correlation bins.
        h = fitsio.FITS(f)
        #nobj += h[1].read_header()['NOBJ']
        we_aux = h[2]['WE'][:]
        wet_aux = we_aux.sum(axis=0)
        rp += h[1]['RP'][:]*wet_aux
        rt += h[1]['RT'][:]*wet_aux
        z  += h[1]['Z'][:]*wet_aux
        nb += h[1]['NB'][:]
        wet += wet_aux
        f_hid = h[2]['HEALPID'][:]

        #Check that the HEALPix pixels are the same.
        if sp.sum(f_hid == hid) == hid.shape[0]:
            da += h[2]['DA'][:] * we_aux
            we += h[2]['WE'][:]
        elif set(f_hid) == set(hid):
            # TODO: Add in check to see if they're the same but just ordered differently.
            raise IOError('Correlations\' pixels are not ordered in the same way!')
        else:
            raise IOError('Correlations do not have the same footprint!')

        h.close()

    # Normalise all variables by the total weights.
    w = we>0
    da[w] /= we[w]
    wt = wet>0
    rp[wt] /= wet[wt]
    rt[wt] /= wet[wt]
    z[wt] /= wet[wt]

    # Update header's nobj:
    #head['NOBJ'] = nobj

    if fout is not None:

        out = fitsio.FITS(fout,'rw',clobber=True)
        out.write([rp,rt,z,nb],names=['RP','RT','Z','NB'],
            comment=['R-parallel','R-transverse','Redshift','Number of pairs'],
            units=['h^-1 Mpc','h^-1 Mpc','',''],
            header=head,extname='ATTRI')

        head2 = [{'name':'HLPXSCHM','value':'RING','comment':'Healpix scheme'}]
        out.write([hid,we,da],names=['HEALPID','WE','DA'],
            comment=['Healpix index', 'Sum of weight', 'Correlation'],
            header=head2,extname='COR')

        out.close()

    return rp,rt,nb,z,wet,da,we,head

for r in v_realisations:
    print('Constructing final random lya_qso_cross for v9.0.{}...'.format(r))
    for zbin in zbins:
        print(' -> looking at zbin {}'.format(zbin))
        #Move the incorrectly labelled correlation.
        command = 'mv {}/analysis/correlation_functions/v9.0.{}/measurements/lya_qso_cross/correlations/xcf_lya_qso_cross_R_{}_{}.fits.gz {}/analysis/correlation_functions/v9.0.{}/measurements/lya_qso_cross/correlations/xcf_lya_qso_cross_R_{}_{}_1.fits.gz'.format(basedir,r,zbin[0],zbin[1],basedir,r,zbin[0],zbin[1])
        subprocess.call(command,shell=True)

        #Combine the two correlations.
        fi = []
        fi += ['{}/analysis/correlation_functions/v9.0.{}/measurements/lya_qso_cross/correlations/xcf_lya_qso_cross_R_{}_{}_1.fits.gz'.format(basedir,r,zbin[0],zbin[1])]
        fi += ['{}/analysis/correlation_functions/v9.0.{}/measurements/lya_qso_cross/correlations/xcf_lya_qso_cross_R_{}_{}_2.fits.gz'.format(basedir,r,zbin[0],zbin[1])]
        fout = '{}/analysis/correlation_functions/v9.0.{}/measurements/lya_qso_cross/correlations/xcf_lya_qso_cross_R_{}_{}.fits.gz'.format(basedir,r,zbin[0],zbin[1])

        rp,rt,nb,z,wet,da,we,head = coadd_correlations(fi,fout=fout)
