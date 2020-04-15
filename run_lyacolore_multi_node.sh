#!/bin/bash -l

################################################################################
## This is a script to run LyaCoLoRe across multiple nodes, useful for reducing
## wall computing time.

################################################################################
## USER DEFINED PARAMS.

# Set the config file that we want to point to. Look at the config files to get
# some more detail on the arguments chosen, and  for instructions on how to find
# out more about the options.
CONFIG_FILE="./input_files/config_files/config_v9.0.ini"

# Set where your CoLoRe output is located, and where you would like your
# LyaCoLoRe output to be located.
COLORE_OUT_LOC="..."
LYACOLORE_OUT_LOC="..."

# Specify the queue to use for this set of realisations, and the number of
# nodes/cores.
QUEUE='debug'
NNODES=32
NCORES=64
TIME="00:10:00" #hh:mm:ss

# Specify the settings for LyaCoLoRe.
RUN_FILE="${LYACOLORE_OUT_LOC}/run_lyacolore.sh"

## END OF USER DEFINED PARAMS.
################################################################################
## Echo the settings and outputs chosen for each realisation.

echo " "
echo "################################################################################"
echo " "
echo "CoLoRe input will be taken from "$COLORE_OUT_LOC
INPUT_FILES=`ls -1 ${COLORE_OUT_LOC}/out_srcs_*.fits`
NFILES=`echo $INPUT_FILES | wc -w`
echo " -> ${NFILES} input files have been found"
echo "Output will written to "$LYACOLORE_OUT_LOC
if [ ! -d $LYACOLORE_OUT_LOC ] ; then
    mkdir -p $LYACOLORE_OUT_LOC
fi
echo " -> Output logs will be saved to "$LYACOLORE_OUT_LOC"/logs"
if [ ! -d $LYACOLORE_OUT_LOC/logs ] ; then
    mkdir -p $LYACOLORE_OUT_LOC/logs
fi

################################################################################
## Make the master file, and the new file structure.
echo " "
echo "Starting LyaCoLoRe..."
echo " "
echo " 1. Make master file"
echo " "
${LYACOLORE_PATH}/scripts/make_master.py -c ${CONFIG_FILE} -i ${COLORE_OUT_LOC} -o ${LYACOLORE_OUT_LOC} --nproc ${NCORES}

################################################################################
## Generate the run files.
echo " "
echo " 2. Make a slurm run file to make transmission skewers"
echo " "
echo "Run file will be written to "$RUN_FILE
cat > $RUN_FILE <<EOF
#!/bin/bash -l

#SBATCH --partition ${QUEUE}
#SBATCH --nodes ${NNODES}
#SBATCH --time ${TIME}
#SBATCH --job-name lyacolore
#SBATCH --error "${LYACOLORE_OUT_LOC}/lyacolore-%j.err"
#SBATCH --output "${LYACOLORE_OUT_LOC}/lyacolore-%j.out"
#SBATCH -C haswell
#SBATCH -A desi

conda activate lyacolore
umask 0002
export OMP_NUM_THREADS=64

PIXDIRS=\`\ls -tr1d ${LYACOLORE_OUT_LOC}/[0-9]*/* | sort -R\`
NPIXELS=\`echo \$PIXDIRS | wc -w\`
PIXDIRS_list=(\$PIXDIRS)

PIXELS=()
for PIXDIR in \$PIXDIRS ; do
    PIX=\${PIXDIR##*/}
    PIXELS=("\${PIXELS[@]}" \$PIX)
done

NPIXELS_PER_NODE=\$(( (\$NPIXELS + $NNODES - 1)/$NNODES ))
START_INDEX=0
STOP_INDEX=\$(( \$NPIXELS_PER_NODE - 1 ))
FINAL_NODE=0

for NODE in \`seq $NNODES\` ; do
    echo "starting node \$NODE"

    NODE_PIXELS=\${PIXELS[@]:\$START_INDEX:\$NPIXELS_PER_NODE}

    echo "looking at pixels: \${NODE_PIXELS}"

    command="${LYACOLORE_PATH}/scripts/make_transmission.py -c ${CONFIG_FILE} -i ${COLORE_OUT_LOC} -o ${LYACOLORE_OUT_LOC} --nproc ${NCORES} --pixels \${NODE_PIXELS}"

    echo \$command
    \$command >& ${LYACOLORE_OUT_LOC}/logs/node-\${NODE}.log &

    if (( \$FINAL_NODE == 1)) ; then
        echo "all pixels allocated, no more nodes needed"
        break
    fi

    START_INDEX=\$(( \$STOP_INDEX + 1 ))
    STOP_INDEX=\$(( \$START_INDEX + \$NPIXELS_PER_NODE - 1))

    if (( \$STOP_INDEX >= (\$NPIXELS - 1) )) ; then
        STOP_INDEX=\$NPIXELS-1
        FINAL_NODE=1
    fi

done
wait
date

EOF

################################################################################
## Send the job to the queue.
echo " "
echo " 3. Send it to the queue"
echo " "
sbatch $RUN_FILE
echo " "
echo "Done!"
echo " "
echo "################################################################################"

# Copy the config file to the output location for clarity.
cp $CONFIG_FILE $LYACOLORE_OUT_LOC

################################################################################
