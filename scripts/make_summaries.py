import numpy as np
from astropy.io import fits

from pyacolore import DLA

parser.add_argument('--base-dir', type = str, default = None, required=True,
                    help = 'base data directory')

parser.add_argument('--nside', type = int, default = 16, required=False,
                    help = 'HEALPix nside for output files (must be 2^n)')

parser.add_argument('--pixels', type = int, default = None, required=False,
                    help = 'which pixel numbers to work on', nargs='*')

base_dir = args.base_dir
N_side = args.nside
pixels = args.pixels
if not pixels:
    pixels = list(range(12*N_side**2))

DLA.make_DLA_master(base_dir,N_side,pixels)
