
################################################################################
# clear out the old example data
rm -r $LYACOLORE_PATH/example_data/lya_skewers/*

# choose which pixels to run on
# note, only pixels xxx to yyy are included in the example data
PIXELS=`seq 0 10`

# specify number of cores to use
NCORES=1

################################################################################

# specify input file information
FILE_FORMAT='colore'
SKEWER_TYPE='gaussian'

# specify process parameters
NSIDE=16
IVAR_CUT=1150.0
CELL_SIZE=0.25
LAMBDA_MIN=3550.0
MIN_CAT_Z=1.8
LYACOLORE_SEED=123
DLA_BIAS=2.0
DLA_BIAS_EVOL='b_const'
DLA_BIAS_METHOD='global'
DOWNSAMPLING=1.0
FOOTPRINT='full_sky'

# specify transmission file properties
TRANS_LMIN=3470.0
TRANS_LMAX=6500.0
TRANS_DL=0.2
TRANS_FORMAT='final'

# specify process flags
MM_FLAGS="--overwrite --pixels $PIXELS"
MT_FLAGS="--add-DLAs --add-RSDs --add-small-scale-fluctuations --overwrite --compress --add-Lyb --add-metals"
MS_FLAGS="--overwrite --compress --compressed-input"

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
TUNING_PATH="${LYACOLORE_PATH}/input_files/tuning_data_with_bias_vel1.3_b1.65_lr1200.fits"

# make master file and new file structure
echo "making master file"
command="${PROCESS_PATH}/make_master.py --in-dir ${INPUT_PATH} --out-dir ${OUTPUT_PATH} --file-format ${FILE_FORMAT} --skewer-type ${SKEWER_TYPE} --nside ${NSIDE} --nproc ${NCORES} --min-cat-z ${MIN_CAT_Z} ${MM_FLAGS} --footprint ${FOOTPRINT} --pixels ${PIXELS}"
$command

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
command="${PROCESS_PATH}/make_transmission.py --in-dir ${INPUT_PATH} --out-dir ${OUTPUT_PATH}  --file-format ${FILE_FORMAT} --skewer-type ${SKEWER_TYPE} ${MT_FLAGS} --pixels ${PIXELS} --tuning-file ${TUNING_PATH} --nside ${NSIDE} --nproc ${NCORES} --IVAR-cut ${IVAR_CUT} --cell-size ${CELL_SIZE} --lambda-min ${LAMBDA_MIN} --seed ${LYACOLORE_SEED} --DLA-bias ${DLA_BIAS} --DLA-bias-evol ${DLA_BIAS_EVOL} --DLA-bias-method ${DLA_BIAS_METHOD} --transmission-lambda-min ${TRANS_LMIN} --transmission-lambda-max ${TRANS_LMAX} --transmission-delta-lambda ${TRANS_DL} --transmission-format ${TRANS_FORMAT}"
$command

echo "producing analysis pixels"
command="${PROCESS_PATH}/make_summaries.py --base-dir ${OUTPUT_PATH} --nproc ${NCORES} --pixels ${PIXELS} --overwrite --picca-N-merge-values 1 10 --compressed-input --compress ${MS_FLAGS}"
#$command
