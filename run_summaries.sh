# specify number of nodes and cores to use
QUEUE='regular'
NCORES=64
TIME="01:30:00" #hh:mm:ss

# specify process parameters
NSIDE=16

# specify N_merge n_values
NMERGE_VALUES='1 10'

# full path to proces_colore executable (parallel version)
PROCESS_PATH="/global/homes/j/jfarr/Projects/LyaCoLoRe/scripts/"

# full path to folder where input will be taken from
BASE_DIR="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v5/v5.0.0/"
echo "input will be taken from "$BASE_DIR
INPUT_FILES=`ls -1 ${BASE_DIR}/*/*/transmission*.fits`

NFILES=`echo $INPUT_FILES | wc -w`
echo "${NFILES} input files have been found"

RUN_FILE="/global/homes/j/jfarr/Projects/LyaCoLoRe/run_files/run_process_colore_summaries.sh"

cat > $RUN_FILE <<EOF
#!/bin/bash -l

#SBATCH --partition ${QUEUE}
#SBATCH --nodes 1
#SBATCH --time ${TIME}
#SBATCH --job-name process_colore_summaries
#SBATCH --error "/global/homes/j/jfarr/Projects/LyaCoLoRe/run_files/process-colore-summaries-%j.err"
#SBATCH --output "/global/homes/j/jfarr/Projects/LyaCoLoRe/run_files/process-colore-summaries-%j.out"
#SBATCH -C haswell
#SBATCH -A desi

umask 0002
export OMP_NUM_THREADS=64

command="srun -N 1 -n 1 -c ${NCORES} ${PROCESS_PATH}/make_summaries.py --base-dir ${BASE_DIR} --nside ${NSIDE} --nproc ${NCORES} --picca-N-merge-values ${NMERGE_VALUES}"

echo \$command
\$command

date

EOF

sbatch $RUN_FILE
