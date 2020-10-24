#!/bin/bash -l

################################################################################
## This is a script to generate one full realisation of LyaCoLoRe's output,
## using one multi-core node.

################################################################################
## USER DEFINED PARAMS.

# Set the config file that we want to point to. Look at the config files to get
# some more detail on the arguments chosen, and  for instructions on how to find
# out more about the options.
CONFIG_FILE="./input_files/config_files/config_v9.0.ini"

# Set where your CoLoRe output is located, and where you would like your
# LyaCoLoRe output to be located.
COLORE_OUT_LOC=...
LYACOLORE_OUT_LOC=...

# Specify number of cores to use.
NCORES=64

## END OF USER DEFINED PARAMS.
################################################################################

echo "Starting LyaCoLoRe..."
echo " "

# Make master file and new file structure
echo " "
echo " 1. Make master file"
echo " "
command="make_master.py -c ${CONFIG_FILE} -i ${COLORE_OUT_LOC} -o ${LYACOLORE_OUT_LOC} --nproc ${NCORES}"
$command

# Make transmission files and other associated skewer files.
echo " "
echo " 2. Make transmission files"
echo " "
command="make_transmission.py -c ${CONFIG_FILE} -i ${COLORE_OUT_LOC} -o ${LYACOLORE_OUT_LOC} --nproc ${NCORES}"
$command

echo " "
echo "Done!"
echo " "

# Copy the config file to the output location for clarify.
cp $CONFIG_FILE $LYACOLORE_OUT_LOC

################################################################################
