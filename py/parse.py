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
                help = 'Input data directory.')
    parser.add('-o','--out-dir', type = str, required=True,
                help = 'Output data directory.')

    ## Input arguments.
    parser.add('--param-file', type = str, required=False,
                help = 'CoLoRe parameter file name.')
    parser.add('--file-format', type = str, required=False,
                choices=['colore'],
                help = 'Format of the input data files.')
    parser.add('--skewer-type', type = str, required=False,
                choices=['gaussian','density'],
                help = 'Quantity stored in input skewers.')

    ## Processing arguments: survey
    parser.add('--min-cat-z', type = float, required=False,
                help = 'Minimum z of objects in catalogue.')
    parser.add('--pixels', type = int, required=False,
                help = 'List of HEALPix pixel numbers to work on.', nargs='*')
    parser.add('--footprint', type = str, required = False,
                choices=['full_sky','desi','desi_pixel','desi_pixel_plus'],
                help = 'Choice of footprint of QSOs on sky.')
    parser.add('--downsampling', type = float, required=False,
                help = 'Proportion by which to downsample the CoLoRe output.')

    ## Processing arguments: skewers
    parser.add('--rest-frame-weights-cut', type = float, required=False,
                help = 'Maximum rest-frame wavelength for pixels to have weight=1 (Angstroms).')
    parser.add('--cell-size', type = float, required=False,
                help = 'Size of the LyaCoLoRe cells to use (Mpc/h).')
    parser.add('--lambda-min', type = float, required=False,
                help = 'Minimum wavelength value in LyaCoLoRe skewers (Angstroms).')
    parser.add('--tuning-file', type = str, required=False,
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

    ## Processing arguments: DLAs
    parser.add('--add-DLAs', action='store_true', required=False,
                help = 'Choice whether to add DLAs to the transmission file or not.')
    parser.add('--DLA-bias', type = float, required=False,
                help = 'Bias of DLAs relative to matter density field (no units).')
    parser.add('--DLA-bias-evol', type = str, required=False,
                choices=['b_const','bD_const'],
                help = 'DLA bias evolution with redshift.')
    parser.add('--DLA-bias-method', type = str, required=False,
                choices=['global','sample'],
                help = 'DLA bias is determined by the global or sample value of sigma_G.')

    ## Processing arguments: misc
    parser.add('--nproc', type = int, required=False,
                help = 'Number of processes to use.')
    parser.add('--nside', type = int, required=False,
                help = 'HEALPix nside for output files (must be 2^n).')
    parser.add('--seed', type = int, required=False,
                help = 'Seed from which to generate random numbers.')

    ## Output arguments
    parser.add('--add-picca-drqs', action='store_true', required=False,
                help = 'Choice whether to save picca format drq files or not.')
    parser.add('--picca-all-absorbers', action='store_true', required=False,
                help = 'Choice whether to add all absorbers to picca output files or not.')
    parser.add('--transmission-only', action='store_true', required=False,
                help = 'Choice whether to only save transmission files or not.')
    parser.add('--transmission-lambda-min', type = float, required=False,
                help = 'Minimum wavelength of skewers in transmission files (Angstroms).')
    parser.add('--transmission-lambda-max', type = float, required=False,
                help = 'Maximum wavelength of skewers in transmission files (Angstroms).')
    parser.add('--transmission-delta-lambda', type = float, required=False,
                help = 'Pixel size of transmission files\' wavelength grid (Angstroms).')
    parser.add('--transmission-format', type = str, required=False,
                choices=['develop','final','single_HDU'],
                help = 'Format of skewers in transmission files.')
    parser.add('--overwrite', action='store_true', required=False,
                help = 'Choice whether to overwrite existing files or not.')
    parser.add('--compress', action='store_true', required=False,
                help = 'Choice whether to compress the output files to .fits.gz or not.')

    args = parser.parse_args()

    return args
