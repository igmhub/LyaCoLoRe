# specify number of nodes and cores to use
QUEUE='regular'
NNODES=48
NCORES=64
TIME="02:00:00" #hh:mm:ss

# specify process parameters
NSIDE=16
IVAR_CUT=1150.0
CELL_SIZE=0.25
LAMBDA_MIN=3550.0
MIN_CAT_Z=1.8

# specify process flags
FLAGS="--add-RSDs --add-DLAs --add-Lyb --add-metals"

# specify details of colore output
COLORE_NGRID=4096
COLORE_NODES=32
R_SMOOTH=2.0

# full path to proces_colore executable (parallel version)
PROCESS_PATH="/global/homes/j/jfarr/Projects/LyaCoLoRe/scripts/"

# full path to folder where input will be taken from
INPUT_PATH="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/output_G_hZsmooth_${COLORE_NGRID}_${COLORE_NODES}_sr${R_SMOOTH}_bm1_biasG18_picos/"
#INPUT_PATH="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/test_input/"
echo "input will be taken from "$INPUT_PATH
INPUT_FILES=`ls -1 ${INPUT_PATH}/out_srcs_s1_*.fits`
NFILES=`echo $INPUT_FILES | wc -w`
echo "${NFILES} input files have been found"

# full path to folder where output will be written
#OUTPUT_PATH="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/process_output_G_hZsmooth_${COLORE_NGRID}_${COLORE_NODES}_sr${R_SMOOTH}_bm1_biasG18_picos_nside${NSIDE}_RSD_lya_lyb/"
OUTPUT_PATH="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v3.3/"
echo "output will written to "$OUTPUT_PATH
if [ ! -d $OUTPUT_PATH ] ; then
    mkdir -p $OUTPUT_PATH
fi
echo "output logs will be saved to "$OUTPUT_PATH"/logs"
if [ ! -d $OUTPUT_PATH/logs ] ; then
    mkdir -p $OUTPUT_PATH/logs
fi

# full path to file with tuning sigma_G data
TUNING_PATH="/global/homes/j/jfarr/Projects/LyaCoLoRe/input_files/tune_small_scale_fluctuations.fits"

# we will create this script
RUN_FILE="/global/homes/j/jfarr/Projects/LyaCoLoRe/run_files/process_colore_output_G_hZsmooth_${COLORE_NGRID}_${COLORE_NODES}_sr${R_SMOOTH}_bm1_biasG18_picos.sh"
echo "run file "$RUN_FILE

# make master file and new file structure
date
echo "making master file"
${PROCESS_PATH}/make_master.py --in-dir ${INPUT_PATH} --out-dir ${OUTPUT_PATH} --nside ${NSIDE} --nproc ${NCORES} --min-cat-z ${MIN_CAT_Z}
wait
date

cat > $RUN_FILE <<EOF
#!/bin/bash -l

#SBATCH --partition ${QUEUE}
#SBATCH --nodes ${NNODES}
#SBATCH --time ${TIME}
#SBATCH --job-name process_colore
#SBATCH --error "/global/homes/j/jfarr/Projects/LyaCoLoRe/run_files/process-colore-%j.err"
#SBATCH --output "/global/homes/j/jfarr/Projects/LyaCoLoRe/run_files/process-colore-%j.out"
#SBATCH -C haswell
#SBATCH -A desi

umask 0002
export OMP_NUM_THREADS=64

PIXDIRS=\`\ls -tr1d ${OUTPUT_PATH}/[0-9]*/*\`
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

    command="srun -N 1 -n 1 -c ${NCORES} ${PROCESS_PATH}/make_transmission.py --in-dir ${INPUT_PATH} --out-dir ${OUTPUT_PATH} --pixels \${NODE_PIXELS} --tuning-file ${TUNING_PATH} --nside ${NSIDE} --nproc ${NCORES} --IVAR-cut ${IVAR_CUT} --cell-size ${CELL_SIZE} --lambda-min ${LAMBDA_MIN} ${FLAGS}"

    echo \$command
    \$command >& ${OUTPUT_PATH}/logs/node-\${NODE}.log &

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

sbatch $RUN_FILE
