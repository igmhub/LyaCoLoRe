import numpy as np
from astropy.io import fits
import fitsio

realisations = range(10)
#corr_types = ['lya_auto','dla_auto','qso_auto','lya_qso_cross','lya_dla_cross']
#namings = ['cf','co','co','xcf','xcf']
corr_types = ['lya_auto','lya_qso_cross','lya_dla_cross']
namings = ['cf','xcf','xcf']

for i,ct in enumerate(corr_types):

    print(ct,namings[i])
    #Set up the data structures to store the information.
    cor_data = []
    cor = {}
    ver = 'v9.0.{}'.format(realisations[0])
    h = fitsio.FITS('/project/projectdirs/desi/users/jfarr/LyaCoLoRe_paper/analysis/correlation_functions/{}/measurements/{}/correlations/{}_{}.fits.gz'.format(ver,ct,namings[i],ct))
    attri = {}
    for k in ['RP','RT','Z','NB']:
        attri[k] = np.zeros(h[1][k][:].shape)

    #Assume that the headers in the different files are all the same (and correct)
    head = h[1].read_header()
    head2 = h[2].read_header()
    h.close()

    #Loop through the files, adding the data to our overarching structures.
    for r in realisations:
        ver = 'v9.0.{}'.format(r)
        print(ver)
        h = fitsio.FITS('/project/projectdirs/desi/users/jfarr/LyaCoLoRe_paper/analysis/correlation_functions/{}/measurements/{}/correlations/{}_{}.fits.gz'.format(ver,ct,namings[i],ct))
        for k in ['RP','RT','Z']:
            attri[k] += h[1][k][:] * h[1]['NB'][:]
        attri['NB'] += h[1]['NB'][:]
        cor_data += [h[2][:]]
        h.close()

    #Ensure that data are correctly normalised.
    for k in ['RP','RT','Z']:
        attri[k] /= attri['NB']

    for k in ['WE','DA']:
        cor[k] = np.concatenate([cd[k] for cd in cor_data],axis=0)

    cor['HEALPID'] = np.concatenate([cd['HEALPID']+realisations[i]*(10**5) for i,cd in enumerate(cor_data)],axis=0)

    fout = '/project/projectdirs/desi/users/jfarr/LyaCoLoRe_paper/analysis/correlation_functions/stack/measurements/{}/correlations/{}_{}_subsamples.fits.gz'.format(ct,namings[i],ct)
    out = fitsio.FITS(fout,'rw',clobber=True)
    out.write([attri[k] for k in ['RP','RT','Z','NB']],names=['RP','RT','Z','NB'],
        comment=['R-parallel','R-transverse','Redshift','Number of pairs'],
        units=['h^-1 Mpc','h^-1 Mpc','',''],
        header=head,extname='ATTRI')

    head2 = [{'name':'HLPXSCHM','value':'RING','comment':'Healpix scheme'}]
    out.write([cor[k] for k in ['HEALPID','WE','DA']],names=['HEALPID','WE','DA'],
        comment=['Healpix index', 'Sum of weight', 'Correlation'],
        header=head2,extname='COR')

    out.close

