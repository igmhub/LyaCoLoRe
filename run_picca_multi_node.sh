#number of pixels
NPIXELS=1000
NPIXPERNODE=80
NCORES=64

# specify number of nodes and cores to use
QUEUE='debug'
NNODES=$((($NPIXELS+$NPIXPERNODE-1) / $NPIXPERNODE))
TIME="00:30:00" #hh:mm:ss

#Set bin properties
RPMIN=0.0
RPMAX=160.0
RTMIN=0.0
RTMAX=160.0
NP=40
NT=40

#Set the fitting bin properties
RMINS='20.0 40.0 60.0'
RMAXS=160.0
AFIXS='fixed free'

ZMIN=2.3
ZEFF=2.4
ZMAX=2.5

#density or flux or gaussian
QUANTITY="flux-rebin-10"

#Quantity code
QC="FF"
CORRTYPE="cf"

#picca nside value
NSIDEPICCA=16

#LyaCoLoRe nside value
NSIDELYACOLORE=16

#smoothing radius
SR=2.0

#set the path to the picca directory
PICCA_PATH="/global/homes/j/jfarr/Programs/picca/"

#decide analysis number
ANALYSIS_ID=35
ANALYSIS_ID=`printf "%03d" ${ANALYSIS_ID}`

#find ID number
PICCA_ID=`ls -1d $PICCA_PATH/picca_analysis_0*/picca* | sort -t / -k 2,2 | tail -1`
PICCA_ID=${PICCA_ID##*_}
PICCA_ID=${PICCA_ID#${PICCA_ID%%[123456789]*}}
PICCA_ID=`printf "%05d" $((PICCA_ID+1))`

#in directory
INDIR="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v7/v7.0.0/"

#Output location
OUTPUT_PATH="${PICCA_PATH}/picca_analysis_${ANALYSIS_ID}/picca_${PICCA_ID}/"
echo "output will written to "$OUTPUT_PATH
if [ ! -d $OUTPUT_PATH ] ; then
    mkdir -p $OUTPUT_PATH
fi
NODE_OUTPUT_PATH="$OUTPUT_PATH/node_cfs/"
mkdir $NODE_OUTPUT_PATH
echo "output logs will be saved to "$OUTPUT_PATH"/logs"
if [ ! -d $OUTPUT_PATH/logs ] ; then
    mkdir -p $OUTPUT_PATH/logs
fi

#Final output file names
OUTPUT_FILE="${OUTPUT_PATH}/cf.fits.gz"
OUTPUT_EXP_FILE="${OUTPUT_PATH}/cf_exp.fits.gz"

# we will create this script
RUN_FILE="${OUTPUT_PATH}/run_picca_cf_${ANALYSIS_ID}_${PICCA_ID}.sh"
echo "run file "${RUN_FILE}

date

cat > $RUN_FILE <<EOF
#!/bin/bash -l

#SBATCH --partition ${QUEUE}
#SBATCH --nodes ${NNODES}
#SBATCH --time ${TIME}
#SBATCH --job-name picca-${CORRTYPE}-${QC}-${NPIXELS}
#SBATCH --error "${OUTPUT_PATH}/picca-cf-%j.err"
#SBATCH --output "${OUTPUT_PATH}/picca-cf-%j.out"
#SBATCH -C haswell
#SBATCH -A desi

umask 0002
export OMP_NUM_THREADS=64

if [ $(( $NPIXPERNODE * $NNODES )) -gt $NPIXELS ]; then
    echo "$(( $NPIXPERNODE * $NNODES )) pixels will be used (rather than ${NPIXELS})."
    echo "This ensures that each of the $NNODES nodes gets the same number of pixels: $NPIXPERNODE."
fi

START_INDEX=0
STOP_INDEX=$(( $NPIXPERNODE - 1 ))
FINAL_NODE=0
N_TOTAL=0
N_REMAINING=${NPIXELS}
i=0

for NODE in \`seq $NNODES\` ; do
    echo "starting node \$NODE"

    N=0
    NODE_FILES=""
    while [ \${N} -lt $NPIXPERNODE ]; do

        F="${INDIR}/\$(( \${i}/100 ))/\${i}/picca-${QUANTITY}-16-\${i}.fits";
        if [ -f "\$F" ]; then 
            NODE_FILES="\$NODE_FILES \$F";
            N=\$(( \$N + 1 ));
            N_TOTAL=\$(( \$N_TOTAL + 1 ));
            N_REMAINING=\$(( \$N_REMAINING - 1 ));
        fi;

        i=\$(( \$i + 1 ));

        if [ \$i -gt $(( 12 * $NSIDELYACOLORE * $NSIDELYACOLORE )) ]; then
            N=$NPIXPERNODE;
            echo "not enough files found for node \$NODE: moving on...";
        fi;

    done

    NODE_OUTPUT_FILE="${NODE_OUTPUT_PATH}/cf_\${START_INDEX}_\${STOP_INDEX}.fits.gz"

    command="srun -N 1 -n 1 -c ${NCORES} ${PICCA_PATH}/bin/picca_cf.py --in-dir \"dummy\" --from-image \${NODE_FILES} --out \${NODE_OUTPUT_FILE} --rp-min ${RPMIN} --rp-max ${RPMAX} --rt-max ${RTMAX} --np ${NP} --nt ${NT} --no-project --nside ${NSIDEPICCA} --nproc 64 --z-cut-min ${ZMIN} --z-cut-max ${ZMAX}"

    echo \$command
    \$command >& ${OUTPUT_PATH}/logs/node-\${NODE}.log &

    START_INDEX=\$(( \$STOP_INDEX + 1 ))
    STOP_INDEX=\$(( \$START_INDEX + $NPIXPERNODE - 1))

    if (( \$FINAL_NODE == 1)) ; then
        echo "all files allocated, no more nodes needed"
        break
    fi

    if (( \${N_TOTAL} >= $(( $NPIXELS - 1)) )) ; then
        FINAL_NODE=1
    fi

done

wait

#Combine all of the individual output files together
/global/homes/j/jfarr/Projects/LyaCoLoRe/picca_analysis/combine_picca_cf_files.py --in-dir ${NODE_OUTPUT_PATH} --out ${OUTPUT_FILE}

#Export
${PICCA_PATH}/bin/picca_export.py --data ${OUTPUT_FILE} --out ${OUTPUT_EXP_FILE}
date

EOF

# we will create a paremter file to store the correlation parameters
PARAMETER_FILE=${OUTPUT_PATH}/parameters.txt
echo "parameterfile "${PARAMETER_FILE}

#make parameter file
cat > $PARAMETER_FILE <<EOF
correl_type = $CORRTYPE
quantity = $QUANTITY
N_side = $NSIDEPICCA
N_pixels = ${NPIXELS}
quantities = ${QC}
rpmin = ${RPMIN}
rpmax = ${RPMAX}
rtmin = ${RTMIN}
rtmax = ${RTMAX}
np = ${NP}
nt = ${NT}
zmin = ${ZMIN}
zmax = ${ZMAX}
EOF

#create config files for doing fits
for RMIN in ${RMINS}; do for RMAX in ${RMAXS}; do for AFIX in ${AFIXS}; do
SUFFIX="_${RMIN%%.*}r_a${AFIX}"
CONFIG_FILE="${OUTPUT_PATH}/config_${CORRTYPE}${SUFFIX}.ini"
source ${PICCA_PATH}/txt_templates/config_${CORRTYPE}.ini $OUTPUT_EXP_FILE $CONFIG_FILE $RPMIN $RPMAX $RTMIN $RTMAX $RMIN $RMAX $AFIX $AFIX
CHI2_FILE="${OUTPUT_PATH}/chi2${SUFFIX}.ini"
source ${PICCA_PATH}/txt_templates/chi2.ini $CHI2_FILE $ZEFF $CORRTYPE $SUFFIX
done; done; done

echo "config and chi2 files written"

# send the job to the queue
sbatch $RUN_FILE
