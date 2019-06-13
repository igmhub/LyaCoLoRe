############################### START OF OPTIONS ###############################

## WHAT TO MEASURE ##
CF=0
XCF=1

## INPUT DATA OPTIONS ##
#number of pixels
NPIXELS=1000
NPIXPERNODE=40
#first quantity to correlate
QUANTITY="flux-rebin-10"
QC="FF"
INDIR="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v7/v7_full_lr1200_tuned_vel1.3/"
#INDIR="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/test_velocity_interpolation_full_no_ssf/"
#second quantity to correlate (if desired)
#QUANTITY2="gaussian"
#INDIR2="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v7/v7.0.0/"
#INDIR2="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/test_velocity_interpolation_full_no_ssf/"
DRQ="${INDIR}/master_picca_RSD.fits"
NSIDELYACOLORE=16

## NERSC OPTIONS ##
# specify number of nodes and cores to use
QUEUE='debug'
NNODES=$((($NPIXELS+$NPIXPERNODE-1) / $NPIXPERNODE))
TIME="00:30:00" #hh:mm:ss
NCORES=64
PICCA_PATH="/global/homes/j/jfarr/Programs/picca/"

## CF OPTIONS ##
#Set bin properties
RPMIN=0.0
RPMAX=160.0
RTMIN=0.0
RTMAX=160.0
NP=40
NT=40
ZMIN=2.4
ZMAX=2.6
ZMINQSO=2.4
ZMAXQSO=2.6
NSIDEPICCA=16

## OUTPUT OPTIONS ##
ANALYSIS_ID=40
CF_OUTPUT_FILENAME="cf.fits.gz"
CF_OUTPUT_EXP_FILENAME="cf_exp.fits.gz"
XCF_OUTPUT_FILENAME="xcf.fits.gz"
XCF_OUTPUT_EXP_FILENAME="xcf_exp.fits.gz"

## FIT OPTIONS ##
#Set the fitting bin properties
RMINS='20.0 40.0 60.0'
RMAXS=160.0
AFIXS='fixed free'

################################ END OF OPTIONS ################################

#format analysis number
ANALYSIS_ID=`printf "%03d" ${ANALYSIS_ID}`

#find ID number
PICCA_ID=`ls -1d $PICCA_PATH/picca_analysis_0*/picca* | sort -t / -k 2,2 | tail -1`
PICCA_ID=${PICCA_ID##*_}
PICCA_ID=${PICCA_ID#${PICCA_ID%%[123456789]*}}
PICCA_ID=`printf "%05d" $((PICCA_ID+1))`

#Output location
OUTPUT_PATH="${PICCA_PATH}/picca_analysis_${ANALYSIS_ID}/picca_${PICCA_ID}/"
echo "output will written to "$OUTPUT_PATH
if [ ! -d $OUTPUT_PATH ] ; then
    mkdir -p $OUTPUT_PATH
fi
NODE_OUTPUT_PATH="$OUTPUT_PATH/node_outputs/"
mkdir $NODE_OUTPUT_PATH
if [ ! -d $OUTPUT_PATH/logs ] ; then
    mkdir -p $OUTPUT_PATH/logs
fi

#make paths to output files
CF_OUTPUT_FILE="${OUTPUT_PATH}/${CF_OUTPUT_FILENAME}"
CF_OUTPUT_EXP_FILE="${OUTPUT_PATH}/${CF_OUTPUT_EXP_FILENAME}"
XCF_OUTPUT_FILE="${OUTPUT_PATH}/${XCF_OUTPUT_FILENAME}"
XCF_OUTPUT_EXP_FILE="${OUTPUT_PATH}/${XCF_OUTPUT_EXP_FILENAME}"

# create scripts
if [ $CF -eq 1 ]; then
CF_RUN_FILE="${OUTPUT_PATH}/run_picca_cf_${ANALYSIS_ID}_${PICCA_ID}.sh"
echo "cf run file "${RUN_FILE}
cat > $CF_RUN_FILE <<EOF
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

INDIR="$INDIR"
INDIR2="$INDIR2"
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
    NODE_FILES2=""
    while [ \${N} -lt $NPIXPERNODE ]; do

        F="${INDIR}/\$(( \${i}/100 ))/\${i}/picca-${QUANTITY}-16-\${i}.fits.gz";

        if [ -z "\${INDIR2}" ]; then

            if [ -f "\$F" ]; then
                NODE_FILES="\$NODE_FILES \$F";
                N=\$(( \$N + 1 ));
                N_TOTAL=\$(( \$N_TOTAL + 1 ));
                N_REMAINING=\$(( \$N_REMAINING - 1 ));
            fi;

        else

            F2="${INDIR2}/\$(( \${i}/100 ))/\${i}/picca-${QUANTITY2}-16-\${i}.fits.gz";

            if [ -f "\$F" ] && [ -f "\$F2" ]; then
                NODE_FILES="\$NODE_FILES \$F";
                NODE_FILES2="\$NODE_FILES2 \$F2";
                N=\$(( \$N + 1 ));
                N_TOTAL=\$(( \$N_TOTAL + 1 ));
                N_REMAINING=\$(( \$N_REMAINING - 1 ));
            fi;

        fi;

        i=\$(( \$i + 1 ));

        if [ \$i -gt $(( 12 * $NSIDELYACOLORE * $NSIDELYACOLORE )) ]; then
            N=$NPIXPERNODE;
            echo "not enough files found for node \$NODE: moving on...";
        fi;

    done

    NODE_OUTPUT_FILE="${NODE_OUTPUT_PATH}/cf_\${START_INDEX}_\${STOP_INDEX}.fits.gz"

    if [ -z "\${INDIR2}" ]; then

        command="srun -N 1 -n 1 -c ${NCORES} ${PICCA_PATH}/bin/picca_cf.py --in-dir \"dummy\" --from-image \${NODE_FILES} --out \${NODE_OUTPUT_FILE} --rp-min ${RPMIN} --rp-max ${RPMAX} --rt-max ${RTMAX} --np ${NP} --nt ${NT} --no-project --nside ${NSIDEPICCA} --nproc 64 --z-cut-min ${ZMIN} --z-cut-max ${ZMAX}";

    else

        command="srun -N 1 -n 1 -c ${NCORES} ${PICCA_PATH}/bin/picca_cf.py --in-dir \"dummy\" --from-image \${NODE_FILES} --in-dir2 \"dummy\" --from-image2 \${NODE_FILES2} --out \${NODE_OUTPUT_FILE} --rp-min ${RPMIN} --rp-max ${RPMAX} --rt-max ${RTMAX} --np ${NP} --nt ${NT} --no-project --nside ${NSIDEPICCA} --nproc 64 --z-cut-min ${ZMIN} --z-cut-max ${ZMAX}";

    fi

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

#Generate the aux files
/global/homes/j/jfarr/Projects/LyaCoLoRe/picca_analysis/make_picca_aux_files.py --base-dir ${OUTPUT_PATH} --corr-type $CORRTYPE --quantity $QUANTITY --npixels ${NPIXELS} --quant-code ${QC} --zmin ${ZMIN} --zmax ${ZMAX} --cf-filename $OUTPUT_FILENAME --cf-exp-filename $OUTPUT_EXP_FILENAME --nside $NSIDEPICCA --rpmin ${RPMIN} --rpmax ${RPMAX} --rtmin ${RTMIN} --rtmax ${RTMAX} --np ${NP} --nt ${NT} --rmin-values ${RMINS} --rmax-values ${RMAXS} --afix-values ${AFIXS}

EOF
fi

if [ $XCF -eq 1 ]; then
XCF_RUN_FILE="${OUTPUT_PATH}/run_picca_xcf_${ANALYSIS_ID}_${PICCA_ID}.sh"
echo "xcf run file "${RUN_FILE}
cat > $XCF_RUN_FILE <<EOF
#!/bin/bash -l

#SBATCH --partition ${QUEUE}
#SBATCH --nodes ${NNODES}
#SBATCH --time ${TIME}
#SBATCH --job-name picca-${CORRTYPE}-${QC}-${NPIXELS}
#SBATCH --error "${OUTPUT_PATH}/picca-xcf-%j.err"
#SBATCH --output "${OUTPUT_PATH}/picca-xcf-%j.out"
#SBATCH -C haswell
#SBATCH -A desi

umask 0002
export OMP_NUM_THREADS=64

if [ $(( $NPIXPERNODE * $NNODES )) -gt $NPIXELS ]; then
    echo "$(( $NPIXPERNODE * $NNODES )) pixels will be used (rather than ${NPIXELS})."
    echo "This ensures that each of the $NNODES nodes gets the same number of pixels: $NPIXPERNODE."
fi

INDIR="$INDIR"
INDIR2="$INDIR2"
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
    NODE_FILES2=""
    while [ \${N} -lt $NPIXPERNODE ]; do

        F="${INDIR}/\$(( \${i}/100 ))/\${i}/picca-${QUANTITY}-16-\${i}.fits.gz";

        if [ -z "\${INDIR2}" ]; then

            if [ -f "\$F" ]; then
                NODE_FILES="\$NODE_FILES \$F";
                N=\$(( \$N + 1 ));
                N_TOTAL=\$(( \$N_TOTAL + 1 ));
                N_REMAINING=\$(( \$N_REMAINING - 1 ));
            fi;

        else

            F2="${INDIR2}/\$(( \${i}/100 ))/\${i}/picca-${QUANTITY2}-16-\${i}.fits.gz";

            if [ -f "\$F" ] && [ -f "\$F2" ]; then
                NODE_FILES="\$NODE_FILES \$F";
                NODE_FILES2="\$NODE_FILES2 \$F2";
                N=\$(( \$N + 1 ));
                N_TOTAL=\$(( \$N_TOTAL + 1 ));
                N_REMAINING=\$(( \$N_REMAINING - 1 ));
            fi;

        fi;

        i=\$(( \$i + 1 ));

        if [ \$i -gt $(( 12 * $NSIDELYACOLORE * $NSIDELYACOLORE )) ]; then
            N=$NPIXPERNODE;
            echo "not enough files found for node \$NODE: moving on...";
        fi;

    done

    NODE_OUTPUT_FILE="${NODE_OUTPUT_PATH}/xcf_\${START_INDEX}_\${STOP_INDEX}.fits.gz"

    command="srun -N 1 -n 1 -c ${NCORES} ${PICCA_PATH}/bin/picca_xcf.py --in-dir \"dummy\" --from-image \${NODE_FILES} --drq ${DRQ} --out \${NODE_OUTPUT_FILE} --rp-min ${RPMIN} --rp-max ${RPMAX} --rt-max ${RTMAX} --np ${NP} --nt ${NT} --no-project --nside ${NSIDEPICCA} --nproc 64 --z-cut-min ${ZMIN} --z-cut-max ${ZMAX} --z-min-obj ${ZMINQSO} --z-max-obj ${ZMAXQSO}";

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

#Generate the aux files
/global/homes/j/jfarr/Projects/LyaCoLoRe/picca_analysis/make_picca_aux_files.py --base-dir ${OUTPUT_PATH} --corr-type $CORRTYPE --quantity $QUANTITY --npixels ${NPIXELS} --quant-code ${QC} --zmin ${ZMIN} --zmax ${ZMAX} --cf-filename $OUTPUT_FILENAME --cf-exp-filename $OUTPUT_EXP_FILENAME --nside $NSIDEPICCA --rpmin ${RPMIN} --rpmax ${RPMAX} --rtmin ${RTMIN} --rtmax ${RTMAX} --np ${NP} --nt ${NT} --rmin-values ${RMINS} --rmax-values ${RMAXS} --afix-values ${AFIXS}

EOF
fi

# send the jobs to the queue
sbatch $CF_RUN_FILE
sbatch $XCF_RUN_FILE
