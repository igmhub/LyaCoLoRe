# LyaCoLoRe
Code development to use CoLoRe simulations for generating simulated Lyman alpha forest spectra

## Install
You only need to include LyaCoLoRe/py in your PYTHONPATH, with something like: 
export PYTHONPATH=$PYTHONPATH:$HOME/Programs/igmhub/LyaCoLoRe/py

## Examples
You can find some examples under example_scripts/. Specifically, the examples available are:

 - plot_delta_picca.ipynb
      To quickly open and view files in the picca input format.

 - plot_raw_colore.ipynb
      To quickly open and view files in the CoLoRe output format.

 - mock_desi.ipynb (requires DESI specific packages)
      To produce mock DESI skewers.

 - read_skewers_lya.ipynb
      To look in detail at files in the CoLoRe output format.

## Main production of DESI mocks

There are two main stages to processing the output files from a CoLoRe simulation, each with a separate script in LyaCoLoRe. These are:

1. Making the master file, and creating the new file structure:

This is carried out by the script 'example_scripts/make_master.py'. In order to run this script, use the following command:

example_scripts/make_master.py --in-dir /path/to/input/directory/ --out-dir /path/to/output/directory/ --nside 16 --nproc 64

The input directory should contain all of the output files from CoLoRe, including the out_params.cfg file. The output directory should be empty. The option "--nside" specifies how fine a HEALPix pixelisation you would like to use (this should be a power of 2). The option "--nproc" specifies the number of processes for the computer to use (a single NERSC compute node can use 64) in order to ensure the script runs quickly.

2. Creating the output files:

This is carried out by the script 'example_scripts/make_transmission.py'. In order to run this script, use the following command:

example_scripts/make_transmission.py --in-dir /path/to/input/directory/ --out-dir /path/to/output/directory/ --nside 16 --nproc 64

The options are as explained in part 1. Other option of interest are:
 - if you are only looking to run on a small number of skewers, the "--pixels" option allows you to specify pixel numbers to work on. For example adding "--pixels 0 1 2 3" would produce output files for pixels 0, 1, 2 and 3 and ignore all other pixels
 - if you would only like to produce transmission files (and not Gaussian or Density files), then the option "--transmission-only" will do this

These two stages can be carried out in parallel using the script 'run_process_colore_multi_node.sh'.

