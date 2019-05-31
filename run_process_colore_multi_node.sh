# specify number of nodes and cores to use
QUEUE='debug'
NNODES=32
NCORES=64
TIME="00:30:00" #hh:mm:ss

# specify process parameters
NSIDE=16
IVAR_CUT=1150.0
CELL_SIZE=0.25
LAMBDA_MIN=3470.0
MIN_CAT_Z=1.8
LYACOLORE_SEED=123
DLA_BIAS=2.0
DLA_BIAS_METHOD='b_const'
DOWNSAMPLING=1.0
VEL_BOOST=1.2
FOOTPRINT='full_sky'
TRANSMISSION_FORMAT='develop'

# specify transmission file wavelength grid
TRANS_LMIN=3470.0
TRANS_LMAX=6500.0
TRANS_DL=0.2

# specify process flags
#MM_FLAGS=""
MM_FLAGS=""
MT_FLAGS="--add-DLAs --add-RSDs --add-QSO-RSDs"

# specify details of colore output
COLORE_NGRID=4096
COLORE_NODES=32
R_SMOOTH=2.0
COLORE_SEED=1003

# full path to proces_colore executable (parallel version)
PROCESS_PATH="/global/homes/j/jfarr/Projects/LyaCoLoRe/scripts/"

# full path to folder where input will be taken from
INPUT_PATH="/project/projectdirs/desi/mocks/lya_forest/develop/london/colore_raw/v5_seed${COLORE_SEED}/"
COLORE_PARAM_PATH="${INPUT_PATH}/param_v5_seed${COLORE_SEED}.cfg"
echo "input will be taken from "$INPUT_PATH
INPUT_FILES=`ls -1 ${INPUT_PATH}/out_srcs_*.fits`

NFILES=`echo $INPUT_FILES | wc -w`
echo "${NFILES} input files have been found"

# code version
V_CODE_MAJ="7"
V_CODE_MIN="3"
V_REALISATION="0"

# full path to folder where output will be written
OUTPUT_PATH="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v${V_CODE_MAJ}/v7_full_no_ssf/"
#OUTPUT_PATH="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v${V_CODE_MAJ}/v${V_CODE_MAJ}.${V_CODE_MIN}.${V_REALISATION}/"
#OUTPUT_PATH="/project/projectdirs/desi/mocks/lya_forest/london/v${V_CODE_MAJ}.${V_CODE_MIN}/v${V_CODE_MAJ}.${V_CODE_MIN}.${V_REALISATION}/"

echo "output will written to "$OUTPUT_PATH
if [ ! -d $OUTPUT_PATH ] ; then
    mkdir -p $OUTPUT_PATH
fi
echo "output logs will be saved to "$OUTPUT_PATH"/logs"
if [ ! -d $OUTPUT_PATH/logs ] ; then
    mkdir -p $OUTPUT_PATH/logs
fi

# full path to file with tuning sigma_G data
TUNING_PATH="/global/homes/j/jfarr/Projects/LyaCoLoRe/input_files/tuning_data_with_bias_vel1.2_b1.65.fits"
#TUNING_PATH="/global/homes/j/jfarr/Projects/LyaCoLoRe/input_files/tuning_data_with_bias_a2.0_b1.65.fits"

# we will create this script
RUN_FILE="/global/homes/j/jfarr/Projects/LyaCoLoRe/run_files/process_colore_v${V_CODE_MAJ}.${V_CODE_MIN}.${V_REALISATION}.sh"
echo "run file "$RUN_FILE

# we create a log of the inputs/parameters used in the run
PARAM_FILE="${OUTPUT_PATH}/input.param"

# make master file and new file structure
date
echo "making master file"
${PROCESS_PATH}/make_master.py --in-dir ${INPUT_PATH} --out-dir ${OUTPUT_PATH} --nside ${NSIDE} --nproc ${NCORES} --min-cat-z ${MIN_CAT_Z} ${MM_FLAGS} --downsampling ${DOWNSAMPLING} --footprint ${FOOTPRINT}
wait
date

cat > $RUN_FILE <<EOF
#!/bin/bash -l

#SBATCH --partition ${QUEUE}
#SBATCH --nodes ${NNODES}
#SBATCH --time ${TIME}
#SBATCH --job-name process_colore
#SBATCH --error "${OUTPUT_PATH}/process-colore-%j.err"
#SBATCH --output "${OUTPUT_PATH}/process-colore-%j.out"
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

    command="srun -N 1 -n 1 -c ${NCORES} ${PROCESS_PATH}/make_transmission.py --in-dir ${INPUT_PATH} --out-dir ${OUTPUT_PATH} ${MT_FLAGS} --pixels \${NODE_PIXELS} --tuning-file ${TUNING_PATH} --nside ${NSIDE} --nproc ${NCORES} --IVAR-cut ${IVAR_CUT} --cell-size ${CELL_SIZE} --lambda-min ${LAMBDA_MIN} --seed ${LYACOLORE_SEED} --DLA-bias ${DLA_BIAS} --DLA-bias-method ${DLA_BIAS_METHOD} --velocity-multiplier ${VEL_BOOST} --transmission-lambda-min ${TRANS_LMIN} --transmission-lambda-max ${TRANS_LMAX} --transmission-delta-lambda ${TRANS_DL} --transmission-format ${TRANSMISSION_FORMAT}"

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

# write the input file
cat > $PARAM_FILE <<EOF
#LyaCoLoRe paths
PROCESS_PATH=${PROCESS_PATH}
INPUT_PATH=${INPUT_PATH}
TUNING_PATH=${TUNING_PATH}
OUTPUT_PATH=${OUTPUT_PATH}

#LyaCoLoRe params
NSIDE=${NSIDE}
IVAR_CUT=${IVAR_CUT}
CELL_SIZE=${CELL_SIZE}
LAMBDA_MIN=${LAMBDA_MIN}
MIN_CAT_Z=${MIN_CAT_Z}
LYACOLORE_SEED=${LYACOLORE_SEED}
DLA_BIAS=${DLA_BIAS}
DLA_BIAS_METHOD=${DLA_BIAS_METHOD}
DOWNSAMPLING=${DOWNSAMPLING}
VEL_BOOST=${VEL_BOOST}
TRANS_LMIN=${TRANS_LMIN}
TRANS_LMAX=${TRANS_LMAX}
TRANS_DL=${TRANS_DL}

#LyaCoLoRe flags
MM_FLAGS=${MM_FLAGS}
MT_FLAGS=${MT_FLAGS}

#CoLoRe params
COLORE_PARAM_PATH=${COLORE_PARAM_PATH}
EOF
cat ${COLORE_PARAM_PATH} >> ${PARAM_FILE}

# copy run file to the output location for record
cp $RUN_FILE $OUTPUT_PATH

# send the job to the queue
sbatch $RUN_FILE
