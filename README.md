# LyaCoLoRe
Code development to use CoLoRe simulations for generating simulated Lyman alpha forest spectra.
It is imortant to notice that the output files of CoLoRe are the input ones for LyaCoLoRe, and LyaCoLoRe's outputs will be the transmission files.

## Install
To install, include LyaCoLoRe/py in your PYTHONPATH, with something like: 
export PYTHONPATH=$PYTHONPATH:$HOME/Programs/igmhub/LyaCoLoRe/py

If you would like to add DLAs using the best avaiable code, you'll need to pip install pyigm. But the code should be able to run without it.

## Examples
You can find some examples under example_scripts/. For instance, you can: 
 - plot_colore_skewer.py : plot a density skewer from CoLoRe
 
 - healpix_quasars.py : plot the angular positions of the quasars in a CoLoRe output, coloured by HEALPix pixel.
 
 - plot_transmission.py : plot the transmitted flux fraction (0 < F < 1) for a processed file.
 
 
## Main production of DESI mocks

There are two main stages to processing the output files from a CoLoRe simulation, each with a separate script in LyaCoLoRe. These are:

1. Making the master file, and creating the new file structure:

This is carried out by the script 'scripts/make_master.py'. In order to run this script, use the following command:

scripts/make_master.py --in-dir example_data/raw_colore_1000/ --out-dir /path/to/output/directory/ --nside 16 --nproc 64

The input directory should contain all of the output files from CoLoRe, including the out_params.cfg file. The output directory should be empty. The option "--nside" specifies how fine a HEALPix pixelisation you would like to use (this should be a power of 2). The option "--nproc" specifies the number of processes for the computer to use (a single NERSC compute node can use 64) in order to ensure the script runs quickly.

Remember that the input directory is a file previously created with CoLoRe. To run this instruction, you should be in the directory where your LyaCoLoRe is cloned.


2. Creating the output files:

This is carried out by the script 'scripts/make_transmission.py'. In order to run this script, use the following command:

scripts/make_transmission.py --in-dir example_data/raw_colore_1000/ --out-dir /path/to/output/directory/ --nside 16 --nproc 64

The options are as explained in part 1. 

Things to keep in mind:
The input directory is a file previously created with CoLoRe. To run this instruction, you should be in the directory where your LyaCoLoRe is cloned. Input and output files are the same as the ones used in stage 1.

Other option of interest are:
 - If you are only looking to run on a small number of skewers, the "--pixels" option allows you to specify pixel numbers to work on. For example adding "--pixels 0 1 2 3" would produce output files for pixels 0, 1, 2 and 3 and ignore all other pixels
 - If you would only like to produce transmission files (and not Gaussian or Density files), then the option "--transmission-only" will do this

These two stages can be carried out in parallel using the script 'run_process_colore_multi_node.sh'.
