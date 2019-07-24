import os
import scipy as sp
import sys
import fitsio
import glob
import healpy
import scipy.interpolate as interpolate
import iminuit
from multiprocessing import Pool
import multiprocessing
import time
import argparse

from picca.data import delta
from lyacolore import utils

################################################################################
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--in-dir', type = str, default = None, required=True,
                    help = 'directory of LyaCoLoRe output files')

parser.add_argument('--out-dir', type = str, default = None, required=True,
                    help = 'directory of output')

parser.add_argument('--in-files', type = str, default = None, required=False,
                    help = 'input files', nargs='*')

parser.add_argument('--nproc', type = int, default = 1, required=False,
                    help = 'number of processes to use')

parser.add_argument('--downsampling', type = float, default = 1.0, required=False,
                    help = 'proportion by which to downsample')

parser.add_argument('--downsampling-seed', type = int, default = 0, required=False,
                    help = 'seed for the random downsampling')

parser.add_argument('--downsample-randoms', action="store_true", default = False, required=False,
                    help = 'whether to downsample the random catalogs or not')

parser.add_argument('--randoms-dir', type = str, default = None, required=False,
                    help = 'directory of randoms')

parser.add_argument('--DLA-z-buffer', type = float, default = 0., required=False,
                    help = 'buffer interval to leave redward of QSOs when making DLA catalogs')

parser.add_argument('--add-Lyb', action="store_true", default = False, required=False,
                    help = 'whether to add Lyb absorption or not')

parser.add_argument('--add-metals', action="store_true", default = False, required=False,
                    help = 'whether to add metal absorption or not')

parser.add_argument('--transmission-lambda-min', type = float, default = 3600., required=False,
                    help = 'minimum wavelength stored in the transmission files')

parser.add_argument('--transmission-lambda-max', type = float, default = 5500., required=False,
                    help = 'maximum wavelength stored in the transmission files')

parser.add_argument('--transmission-lambda-rest-min', type = float, default = 1040., required=False,
                    help = 'minimum wavelength in the rest frame stored in the transmission files')

parser.add_argument('--transmission-lambda-rest-max', type = float, default = 1200., required=False,
                    help = 'maximum wavelength in the rest frame stored in the transmission files')

parser.add_argument('--transmission-delta-lambda', type = float, default = 0.0003, required=False,
                    help = 'pixel size of transmission files wavelength grid')

parser.add_argument('--single-DLA-per-skw', action="store_true", default = False, required=False,
                    help = 'whether to allow at most 1 DLA per skewer or not')

args = parser.parse_args()


# Wavelength grid for output.
lObs_min = args.transmission_lambda_min
lObs_max = args.transmission_lambda_max
lRF_min = args.transmission_lambda_rest_min
lRF_max = args.transmission_lambda_rest_max
dll = args.transmission_delta_lambda

################################################################################

"""
Define the multiprocessing tracking functions
"""

#Define a progress-tracking function.
def log_result(retval):

    results.append(retval)
    N_complete = len(results)
    N_tasks = len(tasks)

    utils.progress_bar(N_complete,N_tasks,start_time)

#Define an error-tracking function.
def log_error(retval):
    print('Error:',retval)

################################################################################
#Make the zcat.

def create_cat(indir,outdir,downsampling,seed=0,DLA_z_buffer=0.,ds_randoms=False,ds_DLA_randoms=False):

    ###
    zmin = 1.70
    nside = 16
    zint = ['0:10']

    ### Make random generator
    state = sp.random.RandomState(seed)

    ### Data
    h = fitsio.FITS(indir+'/master.fits')
    data = {}
    for k in ['RA','DEC']:
        data[k] = h[1][k][:]
    for k in ['THING_ID','PLATE','MJD','FIBERID']:
        data[k] = h[1]['MOCKID'][:]
    data['Z'] = h[1]['Z_QSO_RSD'][:]
    print(data['Z'].min())
    w = data['Z']>zmin
    for k in data.keys():
        data[k] = data[k][w]
    h.close()
    phi = data['RA']*sp.pi/180.
    th = sp.pi/2.-data['DEC']*sp.pi/180.
    pix = healpy.ang2pix(nside,th,phi)
    data['PIX'] = pix
    print('INFO: {} QSO in mocks data'.format(data['RA'].size))

    ### Get reduced data numbers
    original_nbData = data['RA'].shape[0]
    nbData = round(original_nbData * downsampling)

    ### Save data
    assert nbData<=data['RA'].size
    w = state.choice(sp.arange(data['RA'].size), size=nbData, replace=False)
    w_thid = data['THING_ID'][w]
    print('INFO: downsampling to {} QSOs in catalog'.format(nbData))
    out = fitsio.FITS(outdir+'/zcat_{}.fits'.format(downsampling),'rw',clobber=True)
    cols = [ v[w] for k,v in data.items() if k not in ['PIX'] ]
    names = [ k for k in data.keys() if k not in ['PIX'] ]
    out.write(cols,names=names)
    out.close()

    ### DLA data
    h = fitsio.FITS(indir+'/master_DLA.fits')
    data = {}
    for k in ['RA','DEC']:
        data[k] = h[1][k][:]
    for k in ['THING_ID','PLATE','MJD','FIBERID']:
        data[k] = h[1]['MOCKID'][:]
    data['Z'] = h[1]['Z_DLA_RSD'][:]
    # Remove DLAs too close to the QSO
    data['Z_QSO'] = h[1]['Z_QSO_RSD'][:]
    w = abs(data['Z'] - data['Z_QSO'])>DLA_z_buffer
    for k in data.keys():
        data[k] = data[k][w]
    # Remove DLAs below zmin
    w = data['Z']>zmin
    for k in data.keys():
        data[k] = data[k][w]
    h.close()
    phi = data['RA']*sp.pi/180.
    th = sp.pi/2.-data['DEC']*sp.pi/180.
    pix = healpy.ang2pix(nside,th,phi)
    data['PIX'] = pix
    print('INFO: {} DLA in mocks data'.format(data['RA'].size))

    ### Get reduced data numbers
    original_nbData = data['RA'].shape[0]
    nbData = round(original_nbData * downsampling)

    ### Save DLA data
    w_DLA = sp.isin(data['THING_ID'],w_thid)

    if args.single_DLA_per_skw:
        reduced_THING_ID = data['THING_ID'][w_DLA]
        n_id = 1
        current_m = reduced_THING_ID[0]
        ind = 0
        inds = []
        for i,m in enumerate(reduced_THING_ID[1:]):
          i += 1
          if m == current_m:
            n_id += 1
            p = state.uniform()
            if p > 1/n_id:
              ind = i
          else:
            current_m = m
            inds += [ind]
            ind = i
            n_id = 1
        w_DLA = np.array(inds)

    print('INFO: downsampling leaves {} DLAs in catalog'.format(sp.sum(w_DLA)))
    out = fitsio.FITS(outdir+'/zcat_DLA_{}.fits'.format(downsampling),'rw',clobber=True)
    cols = [ v[w_DLA] for k,v in data.items() if k not in ['PIX','Z_QSO'] ]
    names = [ k for k in data.keys() if k not in ['PIX','Z_QSO'] ]
    out.write(cols,names=names)
    out.close()

    if ds_randoms:
        ### Data
        h = fitsio.FITS(indir+'/master_randoms.fits')
        data = {}
        for k in ['RA','DEC']:
            data[k] = h[1][k][:]
        for k in ['THING_ID','PLATE','MJD','FIBERID']:
            data[k] = h[1]['MOCKID'][:]
        data['Z'] = h[1]['Z'][:]
        w = data['Z']>zmin
        for k in data.keys():
            data[k] = data[k][w]
        h.close()
        phi = data['RA']*sp.pi/180.
        th = sp.pi/2.-data['DEC']*sp.pi/180.
        pix = healpy.ang2pix(nside,th,phi)
        data['PIX'] = pix
        print('INFO: {} QSO in randoms'.format(data['RA'].size))

        ### Get reduced data numbers
        original_nbData = data['RA'].shape[0]
        nbData = round(original_nbData * downsampling)

        ### Save data
        assert nbData<=data['RA'].size
        w = state.choice(sp.arange(data['RA'].size), size=nbData, replace=False)
        print('INFO: downsampling to {} QSOs in randoms catalog'.format(nbData))
        out = fitsio.FITS(outdir+'/zcat_{}_randoms.fits'.format(downsampling),'rw',clobber=True)
        cols = [ v[w] for k,v in data.items() if k not in ['PIX'] ]
        names = [ k for k in data.keys() if k not in ['PIX'] ]
        out.write(cols,names=names)
        out.close()

        h = fitsio.FITS(indir+'/master_DLA_randoms.fits')
        data = {}
        for k in ['RA','DEC']:
            data[k] = h[1][k][:]
        for k in ['THING_ID','PLATE','MJD','FIBERID']:
            data[k] = h[1]['MOCKID'][:]
        data['Z'] = h[1]['Z_DLA'][:]
        # Remove DLAs too close to the QSO
        data['Z_QSO'] = h[1]['Z_QSO_RSD'][:]
        w = abs(data['Z'] - data['Z_QSO'])>DLA_z_buffer
        for k in data.keys():
            data[k] = data[k][w]
        # Remove DLAs below zmin
        w = data['Z']>zmin
        for k in data.keys():
            data[k] = data[k][w]
        h.close()
        phi = data['RA']*sp.pi/180.
        th = sp.pi/2.-data['DEC']*sp.pi/180.
        pix = healpy.ang2pix(nside,th,phi)
        data['PIX'] = pix
        print('INFO: {} DLA in randoms'.format(data['RA'].size))

        ### Save DLA data
        w_DLA = sp.isin(data['THING_ID'],w_thid)

        if args.single_DLA_per_skw:
            reduced_THING_ID = data['THING_ID'][w_DLA]
            n_id = 1
            current_m = reduced_THING_ID[0]
            ind = 0
            inds = []
            for i,m in enumerate(reduced_THING_ID[1:]):
              i += 1
              if m == current_m:
                n_id += 1
                p = state.uniform()
                if p > 1/n_id:
                  ind = i
              else:
                current_m = m
                inds += [ind]
                ind = i
                n_id = 1
            w_DLA = np.array(inds)

        print('INFO: downsampling leaves {} DLAs in randoms catalog'.format(sp.sum(w_DLA)))
        out = fitsio.FITS(outdir+'/zcat_DLA_{}_randoms.fits'.format(downsampling),'rw',clobber=True)
        cols = [ v[w_DLA] for k,v in data.items() if k not in ['PIX','Z_QSO'] ]
        names = [ k for k in data.keys() if k not in ['PIX','Z_QSO'] ]
        out.write(cols,names=names)
        out.close()

    return

create_cat(args.in_dir,args.out_dir,args.downsampling,seed=args.downsampling_seed,DLA_z_buffer=args.DLA_z_buffer,ds_randoms=args.downsample_randoms,ds_DLA_randoms=args.downsample_randoms)

################################################################################
#Make the delta files.
### Based on the function desi_convert_transmission_to_delta_files in picca.utils
"""Convert desi transmission files to picca delta files

Args:
    zcat (str): path to the catalog of object to extract the transmission from
    indir (str): path to transmission files directory
    outdir (str): path to write delta files directory
    lObs_min (float) = 3600.: min observed wavelength in Angstrom
    lObs_max (float) = 5500.: max observed wavelength in Angstrom
    lRF_min (float) = 1040.: min Rest Frame wavelength in Angstrom
    lRF_max (float) = 1200.: max Rest Frame wavelength in Angstrom
    dll (float) = 3.e-4: size of the bins in log lambda
    nspec (int) = None: number of spectra, if 'None' use all

Returns:
    None
"""

zcat = args.out_dir+'/zcat_{}.fits'.format(downsampling)
if args.add_Lyb * args.add_metals:
    args.out_dir = args.out_dir+'/deltas_{}_Lyb_metals/'.format(downsampling)
elif args.add_Lyb:
    args.out_dir = args.out_dir+'/deltas_{}_Lyb/'.format(downsampling)
elif args.add_metals:
    args.out_dir = args.out_dir+'/deltas_{}_metals/'.format(downsampling)
else:
    args.out_dir = args.out_dir+'/deltas_{}/'.format(downsampling)
if not os.path.isdir(args.out_dir):
    os.mkdir(args.out_dir)

### Catalog of objects
h = fitsio.FITS(zcat)
key_val = sp.char.strip(sp.array([ h[1].read_header()[k] for k in h[1].read_header().keys()]).astype(str))
if 'TARGETID' in key_val:
    zcat_thid = h[1]['TARGETID'][:]
elif 'THING_ID' in key_val:
    zcat_thid = h[1]['THING_ID'][:]
w = h[1]['Z'][:]>max(0.,lObs_min/lRF_max -1.)
w &= h[1]['Z'][:]<max(0.,lObs_max/lRF_min -1.)
zcat_ra = h[1]['RA'][:][w].astype('float64')*sp.pi/180.
zcat_dec = h[1]['DEC'][:][w].astype('float64')*sp.pi/180.
zcat_thid = zcat_thid[w]
h.close()
print('INFO: Found {} quasars'.format(zcat_ra.size))

### List of transmission files
if (args.in_dir is None and infiles is None) or (args.in_dir is not None and infiles is not None):
    print("ERROR: No transmisson input files or both 'indir' and 'infiles' given")
    sys.exit()
elif args.in_dir is not None:
    fi = glob.glob(args.in_dir+'/*/*/transmission*.fits*')
    fi = sp.sort(sp.array(fi))
    h = fitsio.FITS(fi[0])
    in_nside = h['METADATA'].read_header()['HPXNSIDE']
    nest = h['METADATA'].read_header()['HPXNEST']
    h.close()
    in_pixs = healpy.ang2pix(in_nside, sp.pi/2.-zcat_dec, zcat_ra, nest=nest)
    fi = sp.sort(sp.array(['{}/{}/{}/transmission-{}-{}.fits.gz'.format(args.in_dir,int(f//100),f,in_nside,f) for f in sp.unique(in_pixs)]))
else:
    fi = sp.sort(sp.array(infiles))
print('INFO: Found {} files'.format(fi.size))

### Stack the transmission
lmin = sp.log10(lObs_min)
lmax = sp.log10(lObs_max)
nstack = int((lmax-lmin)/dll)+1

### Read
def get_stack_data(f):
    #Set up variables.
    deltas = {}
    T_stack = sp.zeros(nstack)
    n_stack = sp.zeros(nstack)

    h = fitsio.FITS(f)
    thid = h['METADATA']['MOCKID'][:]
    if sp.in1d(thid,zcat_thid).sum()==0:
        h.close()
    ra = h['METADATA']['RA'][:].astype(sp.float64)*sp.pi/180.
    dec = h['METADATA']['DEC'][:].astype(sp.float64)*sp.pi/180.
    z = h['METADATA']['Z'][:]
    ll = sp.log10(h['WAVELENGTH'].read())
    trans_names = []
    try:
        trans = h['F_LYA'].read()
    except KeyError:
        try:
            trans = h['F'].read()
        except KeyError:
            try:
                trans = h['TRANSMISSION'].read()
            except KeyError:
                raise KeyError('Transmission not found; check file format.')
    if args.add_Lyb:
        try:
            trans_Lyb = h['F_LYB'].read()
            trans *= trans_Lyb
        except KeyError:
            raise KeyError('Lyb transmission not found; only \'final\' format supported currently.')
    if args.add_metals:
        try:
            trans_metals = h['F_METALS'].read()
            trans *= trans_metals
        except KeyError:
            raise KeyError('Metals transmission not found; only \'final\' format supported currently.')
    nObj = z.size
    pixnum = f.split('-')[-1].split('.')[0]

    if trans.shape[0]!=nObj:
        trans = trans.transpose()

    bins = sp.floor((ll-lmin)/dll+0.5).astype(int)
    tll = lmin + bins*dll
    lObs = (10**tll)*sp.ones(nObj)[:,None]
    lRF = (10**tll)/(1.+z[:,None])
    w = sp.zeros_like(trans).astype(int)
    w[ (lObs>=lObs_min) & (lObs<lObs_max) & (lRF>lRF_min) & (lRF<lRF_max) ] = 1
    nbPixel = sp.sum(w,axis=1)
    cut = nbPixel>=50
    cut &= sp.in1d(thid,zcat_thid)
    if cut.sum()==0:
        h.close()

    ra = ra[cut]
    dec = dec[cut]
    z = z[cut]
    thid = thid[cut]
    trans = trans[cut,:]
    w = w[cut,:]
    nObj = z.size
    h.close()

    deltas[pixnum] = []
    for i in range(nObj):
        tll = ll[w[i,:]>0]
        ttrans = trans[i,:][w[i,:]>0]

        bins = sp.floor((tll-lmin)/dll+0.5).astype(int)
        cll = lmin + sp.arange(nstack)*dll
        cfl = sp.bincount(bins,weights=ttrans,minlength=nstack)
        civ = sp.bincount(bins,minlength=nstack).astype(float)

        ww = civ>0.
        if ww.sum()<50: continue
        T_stack += cfl
        n_stack += civ
        cll = cll[ww]
        cfl = cfl[ww]/civ[ww]
        civ = civ[ww]
        deltas[pixnum].append(delta(thid[i],ra[i],dec[i],z[i],thid[i],thid[i],thid[i],cll,civ,None,cfl,1,None,None,None,None,None,None))

    return (n_stack,T_stack,deltas)

tasks = [(f,) for f in fi]

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = nproc)
    results = []
    start_time = time.time()
    for task in tasks:
        pool.apply_async(get_stack_data,task,callback=log_result,error_callback=log_error)
    pool.close()
    pool.join()

### Get stacked transmission
T_stack = sp.zeros(nstack)
n_stack = sp.zeros(nstack)
deltas = {}
for r in results:
    n_stack += r[0]
    T_stack += r[1]
    deltas = {**deltas, **r[2]}

w = n_stack>0.
T_stack[w] /= n_stack[w]

def normalise_deltas(p):
    if len(deltas[p])==0:
        print('No data in {}'.format(p))
    out = fitsio.FITS(args.out_dir+'/delta-{}'.format(p)+'.fits.gz','rw',clobber=True)
    for d in deltas[p]:
        bins = sp.floor((d.ll-lmin)/dll+0.5).astype(int)
        d.de = d.de/T_stack[bins] - 1.
        d.we *= T_stack[bins]**2

        hd = {}
        hd['RA'] = d.ra
        hd['DEC'] = d.dec
        hd['Z'] = d.zqso
        hd['PMF'] = '{}-{}-{}'.format(d.plate,d.mjd,d.fid)
        hd['THING_ID'] = d.thid
        hd['PLATE'] = d.plate
        hd['MJD'] = d.mjd
        hd['FIBERID'] = d.fid
        hd['ORDER'] = d.order

        cols = [d.ll,d.de,d.we,sp.ones(d.ll.size)]
        names = ['LOGLAM','DELTA','WEIGHT','CONT']
        out.write(cols,names=names,header=hd,extname=str(d.thid))
    out.close()
    return

tasks = [(p,) for p in deltas.keys()]

#Run the multiprocessing pool
if __name__ == '__main__':
    pool = Pool(processes = nproc)
    results = []
    start_time = time.time()
    for task in tasks:
        pool.apply_async(normalise_deltas,task,callback=log_result,error_callback=log_error)
    pool.close()
    pool.join()
