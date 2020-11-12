import configargparse
import numpy as np
import os
import sys

def get_args(argv,for_tuning=False):

    try:
        parser = configargparse.ArgParser(default_config_files=[os.environ['LYACOLORE_PATH']+'/input_files/config_files/default.ini'])
    except KeyError:
        raise NameError('Environment variable "LYACOLORE_PATH" has not been defined.')

    ## Config file.
    parser.add('-c','--config', is_config_file=True, required=False,
                help = 'config file')

    ## Required arguments.
    parser.add('-i','--in-dir', type=str, required=True,
                help = 'Input data directory.')
    parser.add('-o','--out-dir', type=str, required=True,
                help = 'Output data directory.')

    ## Input arguments.
    parser.add('--param-file', type=str, required=False,
                help = 'CoLoRe parameter file name.')
    parser.add('--input-filename-prefix', type=str, required=False,
                help = 'Prefix for the CoLoRe files.')
    parser.add('--file-format', type=str, required=False,
                choices=['colore'],
                help = 'Format of the input data files.')
    parser.add('--skewer-type', type=str, required=False,
                choices=['gaussian','density'],
                help = 'Quantity stored in input skewers.')

    ## Processing arguments: survey
    parser.add('--min-cat-z', type=float, required=False,
                help = 'Minimum z of objects in catalogue.')
    parser.add('--pixels', type=int, required=False,
                help = 'List of HEALPix pixel numbers to work on.', nargs='*')
    parser.add('--footprint', type=str, required = False,
                choices=['full_sky','desi','desi_pixel','desi_pixel_plus'],
                help = 'Choice of footprint of QSOs on sky.')
    parser.add('--downsampling', type=float, required=False,
                help = 'Proportion by which to downsample the CoLoRe output.')

    ## Processing arguments: skewers
    parser.add('--rest-frame-weights-cut', type=float, required=False,
                help = 'Maximum rest-frame wavelength for pixels to have weight=1 (Angstroms).')
    parser.add('--cell-size', type=float, required=False,
                help = 'Size of the LyaCoLoRe cells to use (Mpc/h).')
    parser.add('--lambda-min', type=float, required=False,
                help = 'Minimum wavelength value in LyaCoLoRe skewers (Angstroms).')
    parser.add('--tuning-file', type=str, required=False,
                help = 'Name of file from which to load tuning parameters.')
    parser.add('--add-small-scale-fluctuations', action='store_true', required=False,
                help = 'Choice to add small scale fluctuations or not.')
    parser.add('--add-QSO-RSDs', action='store_true', required=False,
                help = 'Choice to add RSDs or not to our QSO redshifts.')
    parser.add('--add-RSDs', action='store_true', required=False,
                help = 'Choice to add RSDs or not to our pixels.')
    parser.add('--include-thermal-effects', action='store_true', required=False,
                help = 'Choice to include thermal effects in RSDs or not')
    parser.add('--add-Lyb', action='store_true', required=False,
                help = 'Choice to add Lyb absorption or not.')
    parser.add('--add-metals', action='store_true', required=False,
                help = 'Choice to add metal absorption or not.')
    parser.add('--metals-selection', type=str, required=False,
                choices = ['standard','full'],
                help = 'Selection of metal absorbers to add to skewers.')
    parser.add('--metals-list', type=str, required=False,
                help = 'List of metal absorbers to add to skewers.', nargs='*')

    ## Processing arguments: DLAs
    parser.add('--add-DLAs', action='store_true', required=False,
                help = 'Choice whether to add DLAs to the transmission file or not.')
    parser.add('--DLA-bias', type=float, required=False,
                help = 'Bias of DLAs relative to matter density field (no units).')
    parser.add('--DLA-bias-evol', type=str, required=False,
                choices=['b_const','bD_const'],
                help = 'DLA bias evolution with redshift.')
    parser.add('--DLA-bias-method', type=str, required=False,
                choices=['global','sample'],
                help = 'DLA bias is determined by the global or sample value of sigma_G.')
    parser.add('--DLA-NHI-min', type=float, required=False,
                help = 'Minimum value of NHI for DLAs.')
    parser.add('--DLA-NHI-max', type=float, required=False,
                help = 'Maximum value of NHI for DLAs.')

    ## Processing arguments: randoms
    parser.add('--rand-factor-QSO', type=float, required=False,
                help = 'Multiplicative factor for number of QSO randoms.')
    parser.add('--rand-method-QSO', type=str, required=False,
                help = 'QSO randoms method.')
    parser.add('--rand-mockid-start', type=int, required=False,
                help = 'Number at which to start random QSO mockids.')
    parser.add('--rand-factor-DLA', type=float, required=False,
                help = 'Multiplicative factor for number of DLA randoms.')
    parser.add('--rand-method-DLA', type=str, required=False,
                help = 'DLA randoms method.')
    parser.add('--add-rand-DLA-NHI', action='store_true', required=False,
                help = 'Include NHI values in random DLAs (slower but more information).')
    parser.add('--rand-dlaid-start', type=int, required=False,
                help = 'Number at which to start random DLAIDs.')

    ## Processing arguments: misc
    parser.add('--nproc', type=int, required=False,
                help = 'Number of processes to use.')
    parser.add('--nside', type=int, required=False,
                help = 'HEALPix nside for output files (must be 2^n).')
    parser.add('--seed', type=int, required=False,
                help = 'Seed from which to generate random numbers.')

    ## Output arguments
    parser.add('--add-picca-drqs', action='store_true', required=False,
                help = 'Choice whether to save picca format drq files or not.')
    parser.add('--picca-all-absorbers', action='store_true', required=False,
                help = 'Choice whether to add all absorbers to picca output files or not.')
    parser.add('--transmission-only', action='store_true', required=False,
                help = 'Choice whether to only save transmission files or not.')
    parser.add('--transmission-lambda-min', type=float, required=False,
                help = 'Minimum wavelength of skewers in transmission files (Angstroms).')
    parser.add('--transmission-lambda-max', type=float, required=False,
                help = 'Maximum wavelength of skewers in transmission files (Angstroms).')
    parser.add('--transmission-delta-lambda', type=float, required=False,
                help = 'Pixel size of transmission files\' wavelength grid (Angstroms).')
    parser.add('--transmission-format', type=str, required=False,
                choices=['develop','final','single_HDU'],
                help = 'Format of skewers in transmission files.')
    parser.add('--overwrite', action='store_true', required=False,
                help = 'Choice whether to overwrite existing files or not.')
    parser.add('--compress', action='store_true', required=False,
                help = 'Choice whether to compress the output files to .fits.gz or not.')

    if not for_tuning:
        args = parser.parse_args(args=argv[1:])
    else:
        args = parser.parse_args(args=argv)

    ## Carry out checks on arguments.
    if np.log2(args.nside)-int(np.log2(args.nside)) != 0:
        raise ValueError('nside must be a power of 2')

    ## If we have density input skewers and want to add DLAs, then raise an
    ## error: this functionality is not yet implemented.
    if (args.skewer_type=='density') & args.add_DLAs:
        raise ValueError('Adding DLAs from density input skewers is not possible yet!')

    # TODO: maybe make this into one argument? just args.metals and then run checks to see what the input is?

    ## If both a selection of metals and a list of metals is supplied then raise
    ## an error: they may be inconsistent.
    if (args.metals_selection is not None) and (args.metals_list is not None):
        raise ValueError('Both a selection of metals and a list of metals have been provided: choose one!')

    # TODO: print to confirm the arguments. e.g. "DLAs will be added"

    return args


def get_tuning_args(argv):

    try:
        parser = configargparse.ArgParser(default_config_files=[os.environ['LYACOLORE_PATH']+'/input_files/tuning_config_files/default.ini'])
    except KeyError:
        raise NameError('Environment variable "LYACOLORE_PATH" has not been defined.')

    ## Config file.
    parser.add('-c','--config', is_config_file=True, required=False,
                help = 'config file')

    ## Required arguments.
    parser.add('--run-config', type=str, required=True,
                help = 'Config file for the LyaCoLoRe run.')
    parser.add('-i','--in-dir', type=str, required=True,
                help = 'Input data directory.')

    ## Tuning initialisation arguments.
    parser.add('--initial-parameter-file', type=str, required=False,
                help = 'File which contains the initial parameter values to start the tuning from.')
    parser.add('--randomise-initial-parameter-values', action='store_true', required=False,
                help = 'Whether to randomise the initial parameter values.')

    ## Implementation arguments.
    parser.add('--n-skewers', type=int, required=False,
                help = 'The number of skewers to use when tuning (no units).')
    parser.add('--z-values', type=float, required=False, nargs='*',
                help = 'Redshift values at which we assess tuning quantities (no units).')
    parser.add('--z-width', type=float, required=False,
                help = 'Width of redshift bins to assess tuning quantities over (no units).')
    parser.add('--k-max', type=float, required=False,
                help = 'Maximum value of k to use when assessing fit to P1D (s km-1).')
    parser.add('--weight-p1d', type=float, required=False,
                help = 'Weighting given to P1D when assessing parameter performance (no units).')
    parser.add('--weight-mean-flux', type=float, required=False,
                help = 'Weighting given to mean flux when assessing parameter performance (no units).')
    parser.add('--weight-bias-delta', type=float, required=False,
                help = 'Weighting given to b_delta when assessing parameter performance (no units).')
    parser.add('--weight-bias-eta', type=float, required=False,
                help = 'Weighting given to b_eta when assessing parameter performance (no units).')
    parser.add('--d-delta', type=float, required=False,
                help = 'Step size to use when calculating b_delta (no units).')
    parser.add('--d-eta', type=float, required=False,
                help = 'Step size to use when calculating b_eta (no units).')
    parser.add('--nproc', type=int, required=False,
                help = 'Number of processes to use.')

    ## Plot arguments.
    parser.add('--k-plot-max', type=float, required=False,
                help = 'Maximum k value to plot in P1D (float, s km-1).')
    parser.add('--plot-dir', type=str, required=False,
                help = 'Directory to which the tuning plots should be saved.')

    ## Output arguments.
    parser.add('--overwrite', action='store_true', required=False,
                help = 'Whether to overwrite an existing output file')

    tuning_args = parser.parse_args()

    ## Print an information line stating that any input directory from the run
    ## config file has been superseded.
    if tuning_args.in_dir is not None:
        print('INFO: Using input directory as specified:\n{}'.format(tuning_args.in_dir))
        print('INFO: Any input directory in the run-config will be ignored.')

    ## Parse the run config file.
    run_args = get_args(
        '-c {} -i {} -o {}'.format(tuning_args.run_config,tuning_args.in_dir,None),
        for_tuning=True
        )

    return tuning_args, run_args

"""
# Data options.
parser.add_argument('--base-dir', type=str, default = None, required=True,
                    help = 'Base directory for the input data')

parser.add_argument('--tuning-file-out', type=str, default = None, required=False,
                    help = 'Out file for the tuning data')

parser.add_argument('--plot-dir-out', type=str, default = None, required=False,
                    help = 'Out directory for the plots')

parser.add_argument('--file-format', type=str, default = 'colore', required=False,
                    choices=['colore'],
                    help = 'input file type')

parser.add_argument('--skewer-type', type=str, default = 'gaussian', required=False,
                    choices=['gaussian','density'],
                    help = 'type of skewer in input file')

# Computational options.
parser.add_argument('--nproc', type=int, default = 1, required=False,
                    help = 'number of processes to use')

# Options for making skewers.
parser.add_argument('--pixels', type=int, default = None, required=False,
                    help = 'Which pixel numbers to use for input files', nargs='*')

parser.add_argument('--N-pixels', type=int, default = None, required=False,
                    help = 'Number of files to use as input')

parser.add_argument('--nside', type=int, default = 16, required=False,
                    help = 'HEALPix nside for output files (must be 2^n)')

parser.add_argument('--min-cat-z', type=float, default = 1.8, required=False,
                    help = 'minimum z of objects in catalog')

parser.add_argument('--seed', type=int, default = 16, required=False,
                    help = 'Random seed to use for the generation of random extra power.')

parser.add_argument('--cell-size', type=float, default = 0.25, required=False,
                    help = 'size in Mpc/h of output cells')

parser.add_argument('--lambda-min', type=float, default = 3550., required=False,
                    help = 'Minimum observed wavelength to use in the tuning process')

parser.add_argument('--lambda-rest-max', type=float, default = 1200., required=False,
                    help = 'Maximum rest-frame wavelength to use in the tuning process')

# Tuning options.
parser.add_argument('--z-values', type=float, default = [2.0,2.4,2.8,3.2], required=False,
                    help = 'which z values to measure at', nargs='*')

parser.add_argument('--z-width', type=float, default = 0.1, required=False,
                    help = 'Width of z bins')

parser.add_argument('--remove-P1D-file', type=str, default = None, required=False,
                    help = 'P1D to remove from the added ssf while tuning')

parser.add_argument('--remove-P1D-measure', action="store_true", default = False, required=False,
                    help = 'Remove the measured P1D from the added ssf while tuning')

parser.add_argument('--fix-all', action="store_true", default = False, required=False,
                    help = 'Fix all parameters.')

# Output options.
parser.add_argument('--k-plot-max', type=float, default = 0.02, required=False,
                    help = 'max value of z to plot')

parser.add_argument('--overwrite', action="store_true", default = False, required=False,
                    help = 'overwrite existing files')

parser.add_argument('--compressed-input', action="store_true", default = False, required=False,
                    help = 'input files in format .fits.gz')

parser.add_argument('--show-plots', action="store_true", default = False, required=False,
                    help = 'show plots')

# TODO: Implement these.
parser.add_argument('--start-from-file', type=str, default = None, required=False,
                    help = 'Tuning data file to use as a start point')

parser.add_argument('--start-from-random', action="store_true", default = False, required=False,
                    help = 'Tuning data file to use as a start point')

parser.add_argument('--nskewers', type=int, default = None, required=False,
                    help = 'number of skewers to process')

parser.add_argument('--downsampling', type=float, default = 1.0, required=False,
                    help = 'fraction by which to subsample the CoLoRe output')
"""
