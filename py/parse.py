import configargparse
import os
import sys

def get_args(argv):

    try:
        parser = configargparse.ArgParser(default_config_files=[os.environ['LYACOLORE_PATH']+'/input_files/config_files/default.ini'])
    except KeyError:
        raise NameError('Environment variable "LYACOLORE_PATH" has not been defined.')

    ## Config file.
    parser.add('-c','--config', is_config_file=True, required=False,
                help = 'config file')

    ## Required arguments.
    parser.add('-i','--in-dir', type = str, required=True,
                help = 'input data directory')
    parser.add('-o','--out-dir', type = str, required=True,
                help = 'output data directory')

    ## Input arguments.
    parser.add('--param-file', type = str, required=False,
                help = 'output parameter file name')

    ## Processing arguments: survey
    parser.add('--min-cat-z', type = float, required=False,
                help = 'minimum z of objects in catalog')
    parser.add('--pixels', type = int, required=False,
                help = 'which pixel numbers to work on', nargs='*')
    parser.add('--footprint', type = str, required = False,
                choices=['full_sky','desi','desi_pixel','desi_pixel_plus'],
                help = 'name of footprint to use')
    parser.add('--downsampling', type = float, required=False,
                help = 'fraction by which to subsample the CoLoRe output')

    ## Processing arguments: skewers
    parser.add('--rest-frame-weights-cut', type = float, required=False,
                help = 'maximum rest frame lambda for IVAR=1 (Å)')
    parser.add('--cell-size', type = float, required=False,
                help = 'size in Mpc/h of output cells')
    parser.add('--lambda-min', type = float, required=False,
                help = 'minimum lambda in picca skewers (Å)')
    parser.add('--tuning-file', type = str, required=False,
                help = 'file name for data about tuning sigma_G/alpha')
    parser.add('--add-small-scale-fluctuations', action='store_true', required=False,
                help = 'add small scale fluctuations to the Gaussian skewers')
    parser.add('--add-RSDs', action='store_true', required=False,
                help = 'add linear RSDs to the transmission file')
    parser.add('--add-Lyb', action='store_true', required=False,
                help = 'add Lyman-beta absorption to TRANSMISSION HDU')
    parser.add('--add-metals', action='store_true', required=False,
                help = 'include metal absorbers in the transmission files')
    parser.add('--include-thermal-effects', action='store_true', required=False,
                help = 'add thermal RSDs to the transmission file')
    parser.add('--add-QSO-RSDs', action='store_true', required=False,
                help = 'add QSO RSDs to the transmission file')

    ## Processing arguments: DLAs
    parser.add('--add-DLAs', action='store_true', required=False,
                help = 'add DLAs to the transmission file')
    parser.add('--DLA-bias', type = float, required=False,
                help = 'bias of DLAs')
    parser.add('--DLA-bias-evol', type = str, required=False,
                choices=['b_const','bD_const'],
                help = 'choose DLA bias evolution with redshift')
    parser.add('--DLA-bias-method', type = str, required=False,
                choices=['global','sample'],
                help = 'choose whether the DLA bias is determined by the global or sample value of sigma_G')

    ## Processing arguments: misc
    parser.add('--nproc', type = int, required=False,
                help = 'number of processes to use')
    parser.add('--nside', type = int, required=False,
                help = 'HEALPix nside for output files (must be 2^n)')
    parser.add('--nskewers', type = int, required=False,
                help = 'number of skewers to process')
    parser.add('--overwrite', action='store_true', required=False,
                help = 'overwrite existing files')
    parser.add('--seed', type = int, required=False,
                help = 'specify seed to generate random numbers')

    ## Output arguments
    parser.add('--add-picca-drqs', action='store_true', required=False,
                help = 'save picca format drq files')
    parser.add('--picca-all-absorbers', action='store_true', required=False,
                help = 'combine all absorbers in the picca tau and flux files')
    parser.add('--transmission-only', action='store_true', required=False,
                help = 'save only the transmission file')
    parser.add('--transmission-lambda-min', type = float, required=False,
                help = 'minimum wavelength stored in the transmission files')
    parser.add('--transmission-lambda-max', type = float, required=False,
                help = 'maximum wavelength stored in the transmission files')
    parser.add('--transmission-delta-lambda', type = float, required=False,
                help = 'pixel size of transmission files wavelength grid')
    parser.add('--transmission-format', type = str, required=False,
                choices=['develop','final','single_HDU'],
                help = 'format of transmission files')
    parser.add('--compress', action='store_true', required=False,
                help = 'compress output files to .fits.gz')

    args = parser.parse_args()

    return args
