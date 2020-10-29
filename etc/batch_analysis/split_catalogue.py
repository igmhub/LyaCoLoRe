import argparse
import os
import numpy as np
from astropy.io import fits

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-c','--catalogue', type = str, default = None, required=True,
                    help = 'catalogue to be split up')

parser.add_argument('-n','--n-split', type = int, default = None, required=True,
                    help = 'number of files to split catalogue into')

parser.add_argument('-o','--out-prefix', type = str, default = None, required=True,
                    help = 'prefix to use for split catalogues')

args = parser.parse_args()

print('Splitting catalogue {} into {} splits'.format(args.catalogue,args.n_split))
catdir = os.path.dirname(args.catalogue)
cat = fits.open(args.catalogue)
nobj_total = len(cat[1].data)
bins = int(args.n_split*np.arange(nobj_total)/nobj_total)

for i in np.arange(args.n_split):
    out = os.path.join([catdir,args.out_prefix+'_{}.fits'.format(i)])
    w = (bins==i)
    print(' -> Writing split {} to file {} with {} objects'.format(i,out,w.sum()))
    cols = fits.ColDefs(cat[1].data[w])
    hdu = fits.BinTableHDU.from_columns(cols,name='CAT')
    hdulist = fits.HDUList([cat[0], hdu)
    hdulist.writeto(out)
    hdulist.close()
