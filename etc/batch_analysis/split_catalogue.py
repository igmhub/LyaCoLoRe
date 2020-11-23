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

parser.add_argument('--seed', type = int, default = 0, required=False,
                    help = 'random seed')

args = parser.parse_args()

print('Splitting catalogue {} into {} splits'.format(args.catalogue,args.n_split))
catdir = os.path.dirname(args.catalogue)

cat = fits.open(args.catalogue)
nobj_total = len(cat[1].data)
gen = np.random.default_rng(seed=args.seed)
bins = (args.n_split*gen.permutation(np.arange(nobj_total))/nobj_total) // 1

for i in np.arange(args.n_split):
    out = os.path.join(catdir, args.out_prefix+'_{}.fits'.format(i))
    w = (bins==i)
    print(' -> Writing split {} to file {} with {} objects'.format(i,out,w.sum()))
    cols = fits.ColDefs(cat[1].data[w])
    hdu = fits.BinTableHDU(cat[1].data[w],name='CAT')
    hdulist = fits.HDUList([cat[0], hdu])
    hdulist.writeto(out)
    hdulist.close()

cat.close()
