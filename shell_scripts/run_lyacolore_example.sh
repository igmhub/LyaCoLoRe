#!/bin/bash -l

################################################################################

## This is a simple script to run LyaCoLoRe on a small dataset of 1000 skewers
## from CoLoRe, included within the repository.

################################################################################
## USER DEFINED PARAMS.

# Set the config file that we want to point to. There are 2 example config
# files available:
#  1. input_files/config_files/example_gaussian.ini
#      - runs off Gaussian output skewers from CoLoRe
#  2. input_files/config_files/example_2lpt.ini
#      - runs off 2LPT output skewers from CoLoRe
# Take a look at the config files to get some more detail on the arguments
# chosen, and for instructions on how to find out more about the options.
CONFIG_FILE=input_files/config_files/example_gaussian.ini

# Specify number of cores to use.
NCORES=1

# TODO: Implement this.
# Set verbosity (1=on, 0=off) for explanation as LyaCoLoRe works.
VERBOSE=1

## END OF USER DEFINED PARAMS.
################################################################################

echo "Starting LyaCoLoRe..."

# Set up output directories if not already there.
for t in gaussian 2lpt; do
  if [ ! -d $LYACOLORE_PATH/example_data/$t/lyacolore_output ] ; then
      mkdir -p $LYACOLORE_PATH/example_data/$t/lyacolore_output
  fi
done

# Make master file and new file structure
echo " "
echo " 1. Make master file"
echo " "
command="make_master.py -c ${CONFIG_FILE} --nproc ${NCORES}"
$command

# Make transmission files and other associated skewer files.
echo " "
echo " 2. Make transmission files"
echo " "
command="make_transmission.py -c ${CONFIG_FILE} --nproc ${NCORES}"
$command

# TODO: Update this.
#echo "producing analysis pixels"
#command="make_summaries.py --base-dir ${OUTPUT_PATH} --nproc ${NCORES} --pixels ${PIXELS} --overwrite --picca-N-merge-values 1 10 --compressed-input --compress ${MS_FLAGS}"
#$command

echo "Done!"
echo " "

################################################################################
