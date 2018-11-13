#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from multiprocessing import Pool
import multiprocessing
import sys
import time
import os
import argparse

from pyacolore import utils, independent, stats, convert, simulation_data, DLA, RSD, tuning

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--base-dir', type = str, default = None, required=True,
                    help = 'directory of LyaCoLoRe output')

parser.add_argument('--nproc', type = int, default = 1, required=False,
                    help = 'number of processes to use')

parser.add_argument('--nside', type = int, default = 16, required=False,
                    help = 'HEALPix nside for output files (must be 2^n)')

parser.add_argument('--pixels', type = int, default = None, required=False,
                    help = 'which pixel numbers to work on', nargs='*')

args = parser.parse_args()

base_dir = args.base_dir
N_processes = args.nproc
N_side = args.nside
pixels = args.pixels
if not pixels:
    pixels = list(range(12*N_side**2))

#Make the DLA master file
DLA.make_DLA_master(base_dir,N_side,pixels)
