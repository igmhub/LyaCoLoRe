# specify number of cores to use
NCORES=64

# specify nside to use
NSIDE=16

# specify remaining parameters
IVAR_CUT=1150.0
CELL_SIZE=0.25
LAMBDA_MIN=3550.0
MIN_CAT_Z=1.8

# specify grid size
NGRID=4096

# specify number of nodes to use
NODES=32

# smoothing radius
R_SMOOTH=2.0

# full path to proces_colore executable (parallel version
PROCESS_PATH="/global/homes/j/jfarr/Projects/LyaCoLoRe/example_scripts/"

# full path to folder where output will be written
INPUT_PATH="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/output_G_hZsmooth_${NGRID}_${NODES}_sr${R_SMOOTH}_bm1_biasG18_picos/"
echo "input will be taken from "$INPUT_PATH

# full path to folder where output will be written
OUTPUT_PATH="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/process_output_G_hZsmooth_${NGRID}_${NODES}_sr${R_SMOOTH}_bm1_biasG18_picos_nside${NSIDE}/"
echo "output will written to "$OUTPUT_PATH
mkdir $OUTPUT_PATH

# full path to file with tuning sigma_G data
TUNING_PATH="/global/homes/j/jfarr/Projects/LyaCoLoRe/input_files/tune_small_scale_fluctuations.fits"

# we will create this script
RUN_FILE="/global/homes/j/jfarr/Projects/LyaCoLoRe/run_files/process_colore_output_G_hZsmooth_${NGRID}_${NODES}_sr${R_SMOOTH}_bm1_biasG18_picos.sh"
echo "run file "$RUN_FILE

cat > $RUN_FILE <<EOF
#!/bin/bash -l

#SBATCH --partition regular
#SBATCH --nodes 1
#SBATCH --time 06:00:00
#SBATCH --job-name process_colore
#SBATCH --error "/global/homes/j/jfarr/Projects/LyaCoLoRe/run_files/process-colore-%j.err"
#SBATCH --output "/global/homes/j/jfarr/Projects/LyaCoLoRe/run_files/process-colore-%j.out"
#SBATCH -C haswell
#SBATCH -A desi

umask 0002
export OMP_NUM_THREADS=64

date
${PROCESS_PATH}/process_colore_multi.py --in-dir ${INPUT_PATH} --out-dir ${OUTPUT_PATH} --tuning-file ${TUNING_PATH} --nside ${NSIDE} --nproc ${NCORES} --IVAR-cut ${IVAR_CUT} --cell-size ${CELL_SIZE} --lambda-min ${LAMBDA_MIN} --min-cat-z ${MIN_CAT_Z} --add-DLAs
date

EOF

sbatch $RUN_FILE

