# specify number of nodes and cores to use
QUEUE='debug'
NNODES=32
NCORES=64
TIME="00:30:00" #hh:mm:ss

# specify process parameters
NSIDE=16

# specify N_merge n_values
NMERGE_VALUES='1 10'

# full path to proces_colore executable (parallel version)
PROCESS_PATH="/global/homes/j/jfarr/Projects/LyaCoLoRe/scripts/"

# full path to folder where input will be taken from
BASE_DIR="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v5/v5.0.0/"
echo "input will be taken from "$INPUT_PATH
INPUT_FILES=`ls -1 ${INPUT_PATH}/*/*/transmission*.fits`

NFILES=`echo $INPUT_FILES | wc -w`
echo "${NFILES} input files have been found"

cat > $RUN_FILE <<EOF
#!/bin/bash -l

#SBATCH --partition ${QUEUE}
#SBATCH --nodes ${NNODES}
#SBATCH --time ${TIME}
#SBATCH --job-name process_colore_summaries
#SBATCH --error "/global/homes/j/jfarr/Projects/LyaCoLoRe/run_files/process-colore-summaries-%j.err"
#SBATCH --output "/global/homes/j/jfarr/Projects/LyaCoLoRe/run_files/process-colore-summaries-%j.out"
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

    command="srun -N 1 -n 1 -c ${NCORES} ${PROCESS_PATH}/make_summaries.py --base-dir ${BASE_DIR} --pixels \${NODE_PIXELS} --nside ${NSIDE} --nproc ${NCORES} --picca-N-merge-values ${NMERGE_VALUES}"

    echo \$command
    \$command

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
