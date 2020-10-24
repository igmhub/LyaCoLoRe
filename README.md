[![Build Status](https://travis-ci.org/igmhub/LyaCoLoRe.svg?branch=tidy_inputs)](https://travis-ci.org/igmhub/LyaCoLoRe)

# LyaCoLoRe
Code development to use CoLoRe simulations for generating simulated Lyman alpha forest spectra.
LyaCoLoRe takes the output files from CoLoRe as an input, carries out several stages of processing, and produces realistic skewers of transmitted flux fraction as an output.
Also in the repository are tools to tune the parameters within LyaCoLoRe's transformation, and to measure the 1D power spectrum of output skewers quickly.

## Install
To install, run:
```bash
python setup.py install --user
```

You then need to include `LyaCoLoRe/py` in your `PYTHONPATH`, by adding:
```bash
export PYTHONPATH=$PYTHONPATH:<path to LyaCoLoRe>/py
```
to your .bashrc file. Then, you should add the path to LyaCoLoRe to your .bashrc file as follows:
```bash
export LYACOLORE_PATH=<path to LyaCoLoRe>
```

If you would like to add DLAs using the best avaiable code, you'll need to install pyigm, the instructions for which are at https://github.com/pyigm/pyigm/blob/master/docs/install.rst. The code will run without it, but the DLA distributions will be more basic.

Other requirements can be installed via:
```bash
pip install requirements.txt
```

## Examples
The simplest way to run LyaCoLoRe is using the bash script `run_lyacolore_example.sh`. This uses a small sample of CoLoRe output (stored within `example_data`), and transformation parameters stored in `input_files/config_files/example.ini`.

You can find some other examples under `example_scripts/`. For instance, you can:

*   `plot_colore_skewer.py`: plot a density skewer from CoLoRe

*   `healpix_quasars.py`: plot the angular positions of the quasars in a CoLoRe output,
    coloured by HEALPix pixel.

*   `plot_transmission.py`: plot the transmitted flux fraction (0 < F < 1) for a processed file.


## Main production of DESI mocks

There are two main stages to processing the output files from a CoLoRe simulation, each with a separate script in LyaCoLoRe. These can be executed with the script `run_lyacolore.sh` (see below for notes on parallisation).

The two main stages of LyaCoLoRe are:

1. Making the master file, and creating the new file structure:

This is carried out by the script `scripts/make_master.py`. In order to run this script, use the following command:
```bash
scripts/make_master.py -c <path to config file>
```

The input directory should contain all of the output files from CoLoRe, including
the ourput parameter file (named `out_params.cfg` normally). The output directory should be empty.

Argument choices can be made via a config file, and there is a template file stored in `input_files/config_files/template.ini` which describes all of the options and gives the data type and units for each one. Default values for arguments (where applicable) are stored in `input_files/config_files/default.ini`.

Further, arguments can also be set via the command line (e.g. `--nside 16`). These will override any arguments set via the config file.

2. Creating the output files:

This is carried out by the script `scripts/make_transmission.py`.
In order to run this script, use the following command:
```bash
scripts/make_transmission.py -c <path to config file>
```

The same config file is used here as above, to ensure consistency.

## Useful options

Other option of interest are:

*   If you are only looking to run on a small number of skewers, the `--pixels` option
allows you to specify pixel numbers to work on. For example adding `--pixels 0 1 2 3`
would produce output files for HEALPix pixels 0, 1, 2 and 3 (using the nested naming convention) and ignore all other pixels

*   If you would only like to produce transmission files (and not Gaussian or Density files), then the option `--transmission-only` will do this

*   You probably want to add RSD to the flux skewers. If so, you'll need to add the
flag `--add-RSDs`.

*   The flags `--add-Lyb` and `--add-metals` will add these to the files.

## Parallelisation

The entirety of the LyaCoLoRe process is broken up into per-HEALPix pixel actions. As such it is (almost) entirely embarassingly parallelised. Normally, this uses `multiprocessing` for parallelising across many cores on a single machine.

If you would like to parallelise across multiple machines/nodes, then the script `run_lyacolore_multi_node.sh` can be used to do so. This runs LyaCoLoRe completely separately for groups of HEALPix pixels, each on a separate node. There is currently no communication between nodes, and is thus very slightly different to running all pixels on 1 node (some quantities are averaged across all pixels while computing the transmission files).
