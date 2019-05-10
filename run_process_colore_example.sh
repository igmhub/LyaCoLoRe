# path to LyaCoLoRe
LYACOLORE_PATH="/Users/James/Projects/LyaCoLoRe/"

# specify number of cores to use
NCORES=2

################################################################################

# specify process parameters
NSIDE=16
IVAR_CUT=1150.0
CELL_SIZE=0.25
LAMBDA_MIN=3550.0
MIN_CAT_Z=1.8
LYACOLORE_SEED=123
DLA_BIAS=2.0
DLA_BIAS_METHOD='b_const'
DOWNSAMPLING=0.5
VEL_BOOST=1.2

# specify transmission file wavelength grid
TRANS_LMIN=3470.0
TRANS_LMAX=6500.0
TRANS_DL=0.2

# specify process flags
MM_FLAGS="--overwrite"
MT_FLAGS="--add-DLAs --add-RSDs --add-QSO-RSDs --add-small-scale-fluctuations --add-Lyb --overwrite"

# specify details of colore output
COLORE_NGRID=4096
COLORE_NODES=32
R_SMOOTH=2.0
COLORE_SEED=1003

# full path to proces_colore executable (parallel version)
PROCESS_PATH="${LYACOLORE_PATH}/scripts/"

# full path to folder where input will be taken from
INPUT_PATH="${LYACOLORE_PATH}/example_data/raw_colore_1000/"
echo "input will be taken from "$INPUT_PATH
INPUT_FILES=`ls -1 ${INPUT_PATH}/out_srcs_*.fits`
NFILES=`echo $INPUT_FILES | wc -w`
echo "${NFILES} input files have been found"

# code version
V_CODE_MAJ="7"
V_CODE_MIN="0"
V_REALISATION="0"

# full path to folder where output will be written
OUTPUT_PATH="${LYACOLORE_PATH}/example_data/lya_skewers/"
echo "output will written to "$OUTPUT_PATH
if [ ! -d $OUTPUT_PATH ] ; then
    mkdir -p $OUTPUT_PATH
fi
echo "output logs will be saved to "$OUTPUT_PATH"/logs"
if [ ! -d $OUTPUT_PATH/logs ] ; then
    mkdir -p $OUTPUT_PATH/logs
fi

# full path to file with tuning sigma_G data
TUNING_PATH="${LYACOLORE_PATH}/input_files/tuning_data_v7.0.0.fits"

# make master file and new file structure
echo "making master file"
command="${PROCESS_PATH}/make_master.py --in-dir ${INPUT_PATH} --out-dir ${OUTPUT_PATH} --nside ${NSIDE} --nproc ${NCORES} --min-cat-z ${MIN_CAT_Z} ${MM_FLAGS}"
$command #|& tee ${OUTPUT_PATH}/logs/node-1.log

PIXDIRS=`ls -tr1d ${OUTPUT_PATH}/[0-9]*/*`
NPIXELS=`echo $PIXDIRS | wc -w`
PIXDIRS_list=($PIXDIRS)

PIXELS=()
for PIXDIR in $PIXDIRS ; do
    PIX=${PIXDIR##*/}
    PIXELS=("${PIXELS[@]}" $PIX)
done

PIXELS=${PIXELS[@]:0:$NPIXELS}

echo "looking at pixels: ${NODE_PIXELS}"
command="${PROCESS_PATH}/make_transmission.py --in-dir ${INPUT_PATH} --out-dir ${OUTPUT_PATH} ${MT_FLAGS} --pixels ${PIXELS} --tuning-file ${TUNING_PATH} --nside ${NSIDE} --nproc ${NCORES} --IVAR-cut ${IVAR_CUT} --cell-size ${CELL_SIZE} --lambda-min ${LAMBDA_MIN} --seed ${LYACOLORE_SEED} --DLA-bias ${DLA_BIAS} --DLA-bias-method ${DLA_BIAS_METHOD} --velocity-multiplier ${VEL_BOOST} --transmission-lambda-min ${TRANS_LMIN} --transmission-lambda-max ${TRANS_LMAX} --transmission-delta-lambda ${TRANS_DL}"
$command #|& tee -a ${OUTPUT_PATH}/logs/node-1.log

echo "producing analysis pixels"
command="${PROCESS_PATH}/make_summaries.py --base-dir ${OUTPUT_PATH} --nproc ${NCORES} --pixels ${PIXELS} --overwrite --picca-N-merge-values 1 10"
$command #|& tee -a ${OUTPUT_PATH}/logs/node-1.log
