import numpy as np
from astropy.io import fits

realisations = range(10)
corr_types = ['lya_auto']

for ct in corr_types:

    #Set up the data structures to store the information.
    cor_data = []
    cor = {}
    h = fits.open('/project/projectdirs/desi/users/jfarr/LyaCoLoRe_paper/analysis/correlation_functions/{}/measurements/{}/correlations/cf_{}.fits.gz'.format(vers[0],ct,ct))
    attri = {}
    for k in ['RP','RT','Z','NB']
        data[k] = np.zeros(h[1].data[k].shape)
    attri_header = fits.Header()
    for k in ['RPMIN','RPMAX','RTMAX','NP','NT','ZCUTMIN','ZCUTMAX','NSIDE']:
        cor_header[k] = h[1].header[k]
    cor_header = fits.Header()
    for k in ['HLPXSCHM']:
        cor_header[k] = h[2].header[k]
    h.close()

    #Loop through the files, adding the data to our overarching structures.
    for r in realisations:
        ver = 'v9.0.{}'.format(r)
        h = fits.open('/project/projectdirs/desi/users/jfarr/LyaCoLoRe_paper/analysis/correlation_functions/{}/measurements/{}/correlations/cf_{}.fits.gz'.format(ver,ct,ct))
        for k in ['RP','RT','Z']:
            attri[k] += h[1].data[k]*h[1].data['NB']
        attri['NB'] += h[1].data['NB']
        cor += [h[2].data]
        h.close()

    #Ensure that data are correctly normalised.
    for k in ['RP','RT','Z']:
        attri[k] /= attri['NB']

    for k in ['WE','DA']:
        cor[k] = np.concatenate([cd[k] for cd in cor_datas])

    cor['HEALPID'] = np.concatenate([cd['HEALPID']+realisations[i]*(10**5) for i,cd in enumerate(cor_data)])

    attri_table = np.array([attri[k]] for k in ['RP','RT','Z','NB'],dtype = ['RP','RT','Z','NB'])
    cor_table = np.array([attri[k]] for k in ['WE','DA','HEALPID'],dtype = ['WE','DA','HEALPID'])

    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)
    hdu_attri = fits.BinTableHDU(attri_table,header=attri_header,name='ATTRI')
    hdu_cor = fits.BinTableHDU(cor_table,header=cor_header,name='COR')
    hdulist = fits.HDUList([prihdu,hdu_attri,hdu_cor])
    hdulist.writeto('/project/projectdirs/desi/users/jfarr/LyaCoLoRe_paper/analysis/correlation_functions/stack/measurements/{}/correlations/cf_{}_subsamples.fits.gz'.format(ver,ct,ct))
    hdulist.close()
