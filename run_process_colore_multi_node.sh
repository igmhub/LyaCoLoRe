# specify number of nodes and cores to use
NNODES=32
NCORES=64

# specify process parameters
NSIDE=16
IVAR_CUT=1150.0
CELL_SIZE=0.25
LAMBDA_MIN=3550.0
MIN_CAT_Z=1.8

# specify process flags
FLAGS="--add-RSDs --add-DLAs"

# specify details of colore output
COLORE_NGRID=4096
COLORE_NODES=32
R_SMOOTH=2.0

# full path to proces_colore executable (parallel version)
PROCESS_PATH="/global/homes/j/jfarr/Projects/LyaCoLoRe/example_scripts/"

# full path to folder where input will be taken from
INPUT_PATH="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/output_G_hZsmooth_${COLORE_NGRID}_${COLORE_NODES}_sr${R_SMOOTH}_bm1_biasG18_picos/"
echo "input will be taken from "$INPUT_PATH
INPUT_FILES=`ls -1 ${INPUT_PATH}/out_srcs_*.fits`
NFILES=`echo $files | wc -w`
echo "${NFILES} input files have been found"

# full path to folder where output will be written
OUTPUT_PATH="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/process_output_G_hZsmooth_${COLORE_NGRID}_${COLORE_NODES}_sr${R_SMOOTH}_bm1_biasG18_picos_nside${NSIDE}/"
echo "output will written to "$OUTPUT_PATH
if [ ! -d $OUTPUT_PATH ] ; then
    mkdir -p $OUTPUT_PATH
fi
echo "output logs will be saved to "$OUTPUT_PATH"/logs"
if [ ! -d $outdir/logs ] ; then
    mkdir -p $outdir/logs
fi

# full path to file with tuning sigma_G data
TUNING_PATH="/global/homes/j/jfarr/Projects/LyaCoLoRe/input_files/tune_small_scale_fluctuations.fits"

# we will create this script
RUN_FILE="/global/homes/j/jfarr/Projects/LyaCoLoRe/run_files/process_colore_output_G_hZsmooth_${COLORE_NGRID}_${COLORE_NODES}_sr${R_SMOOTH}_bm1_biasG18_picos.sh"
echo "run file "$RUN_FILE

cat > $RUN_FILE <<EOF
#!/bin/bash -l

#SBATCH --partition regular
#SBATCH --nodes ${NFILES}
#SBATCH --time 00:30:00
#SBATCH --job-name process_colore
#SBATCH --error "/global/homes/j/jfarr/Projects/LyaCoLoRe/run_files/process-colore-%j.err"
#SBATCH --output "/global/homes/j/jfarr/Projects/LyaCoLoRe/run_files/process-colore-%j.out"
#SBATCH -C haswell
#SBATCH -A desi

umask 0002
export OMP_NUM_THREADS=64

date

echo "making master file"
${PROCESS_PATH}/make_master.py --in-dir ${INPUT_PATH} --out-dir ${OUTPUT_PATH} --nside ${NSIDE} --nproc ${NCORES} --min-cat-z ${MIN_CAT_Z}

date

NPIXELS=$(( 12*$NSIDE*$NSIDE ))
NPIXELS_PER_NODE=$(( $NPIXELS/$NNODES + 1 ))
PIXEL_START=1
PIXEL_STOP=$NPIXELS_PER_NODE

for NODE in `seq $NNODES` ; do
    echo "starting node $NODE"

    # list of files to run
    if (( $NODE == $NNODES )) ; then
        last=""
    fi
    echo ${first}-${last}
    PIXEL_START=$(( $PIXEL_START + $NPIXELS_PER_NODE ))
    PIXEL_STOP=$(( $PIXEL_STOP + $NPIXELS_PER_NODE ))

    #New command for running processing, needs to be updated
    #Takes a master file as input (default location already in output directory)
    #Also takes a list of pixels to deal with
    #Then opens master, finds out where the data for its pixels is stored, and processes on a pixel by pixel basis
    #May need to think about how the input files are accessed

    command="srun -N 1 -n 1 -c ${NCORES} ${PROCESS_PATH}/make_transmission.py --in-dir ${INPUT_PATH} --out-dir ${OUTPUT_PATH} --pix-start ${PIXEL_START} --pix-stop ${PIXEL_STOP} --tuning-file ${TUNING_PATH} --nside ${NSIDE} --nproc ${NCORES} --IVAR-cut ${IVAR_CUT} --cell-size ${CELL_SIZE} --lambda-min ${LAMBDA_MIN} ${FLAGS}"

    echo $command
    $command >& $outdir/logs/node-$node.log &

date

EOF

sbatch $RUN_FILE
