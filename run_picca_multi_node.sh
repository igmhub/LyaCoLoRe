############################### START OF OPTIONS ###############################

## WHAT TO MEASURE ##
CF=0
XCF=1

## INPUT DATA OPTIONS ##
#number of pixels
NPIXELS=1000
NPIXPERNODE=50
#first quantity to correlate
QUANTITY="flux-rebin-10"
INDIR="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v7/v7_full_lr1200_tuned/"
#INDIR="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/test_velocity_interpolation_full_no_ssf/"
#second quantity to correlate (if desired)
#QUANTITY2="gaussian"
#INDIR2="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/v7/v7.0.0/"
#INDIR2="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/test_velocity_interpolation_full_no_ssf/"
DRQ="${INDIR}/master_picca_RSD.fits"

## OTHER PICCA OPTIONS ##
XCFSHUFFLESEED=1
NSIDELYACOLORE=16

## NERSC OPTIONS ##
# specify number of nodes and cores to use
QUEUE='debug'
NNODES=$((($NPIXELS+$NPIXPERNODE-1) / $NPIXPERNODE))
TIME="00:30:00" #hh:mm:ss
NCORES=64
PICCA_PATH="/global/homes/j/jfarr/Programs/picca/"

## CF OPTIONS ##
CF_QC="Fq"
#Set bin properties
CF_RPMIN=0.0
CF_RPMAX=200.0
CF_RTMIN=0.0
CF_RTMAX=200.0
CF_NP=50
CF_NT=50
CF_ZMIN=2.4
CF_ZMAX=2.6
CF_NSIDEPICCA=16

## XCF OPTIONS ##
XCF_QC="Fq"
#Set bin properties
XCF_RPMIN=0.0
XCF_RPMAX=200.0
XCF_RTMIN=-200.0
XCF_RTMAX=200.0
XCF_NP=50
XCF_NT=100
XCF_ZMIN=2.4
XCF_ZMAX=2.6
XCF_ZMINQSO=1.8
XCF_ZMAXQSO=4.0
XCF_NSIDEPICCA=16

## OUTPUT OPTIONS ##
ANALYSIS_ID=41
CF_OUTPUT_FILENAME="cf.fits.gz"
CF_OUTPUT_EXP_FILENAME="cf_exp.fits.gz"
XCF_OUTPUT_FILENAME="xcf.fits.gz"
if [ ! -z $XCFSHUFFLESEED ]; then XCF_SHUFFLE_OUTPUT_FILENAME="xcf_shuffle.fits.gz"; fi
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
PICCA_DIR_EXAMPLE=`ls -1d $PICCA_PATH/picca_analysis_0*/picca* | tail -1`
NSLASH=`echo $PICCA_DIR_EXAMPLE | grep -o "/" | wc -l`
SORT_FIELD=$(( $NSLASH + 1 ))
PICCA_ID=`ls -1d $PICCA_PATH/picca_analysis_0*/picca* | sort -t / -k $SORT_FIELD,$SORT_FIELD | tail -1`
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
if [ ! -z ${XCFSHUFFLESEED} ]; then XCF_SHUFFLE_OUTPUT_FILE="${OUTPUT_PATH}/${XCF_SHUFFLE_OUTPUT_FILENAME}"; fi
XCF_OUTPUT_EXP_FILE="${OUTPUT_PATH}/${XCF_OUTPUT_EXP_FILENAME}"

# create scripts
if [ $CF -eq 1 ]; then
CF_RUN_FILE="${OUTPUT_PATH}/run_picca_cf_${ANALYSIS_ID}_${PICCA_ID}.sh"
echo "cf run file "${CF_RUN_FILE}
cat > $CF_RUN_FILE <<EOF
#!/bin/bash -l

#SBATCH --partition ${QUEUE}
#SBATCH --nodes ${NNODES}
#SBATCH --time ${TIME}
#SBATCH --job-name picca-cf-${QC}-${NPIXELS}
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

        command="srun -N 1 -n 1 -c ${NCORES} ${PICCA_PATH}/bin/picca_cf.py --in-dir \"dummy\" --from-image \${NODE_FILES} --out \${NODE_OUTPUT_FILE} --rp-min ${CF_RPMIN} --rp-max ${CF_RPMAX} --rt-max ${CF_RTMAX} --np ${CF_NP} --nt ${CF_NT} --no-project --nside ${CF_NSIDEPICCA} --nproc 64 --z-cut-min ${CF_ZMIN} --z-cut-max ${CF_ZMAX}";

    else

        command="srun -N 1 -n 1 -c ${NCORES} ${PICCA_PATH}/bin/picca_cf.py --in-dir \"dummy\" --from-image \${NODE_FILES} --in-dir2 \"dummy\" --from-image2 \${NODE_FILES2} --out \${NODE_OUTPUT_FILE} --rp-min ${CF_RPMIN} --rp-max ${CF_RPMAX} --rt-max ${CF_RTMAX} --np ${CF_NP} --nt ${CF_NT} --no-project --nside ${CF_NSIDEPICCA} --nproc 64 --z-cut-min ${CF_ZMIN} --z-cut-max ${CF_ZMAX}";

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
/global/homes/j/jfarr/Projects/LyaCoLoRe/picca_analysis/combine_picca_cf_files.py --in-dir ${NODE_OUTPUT_PATH} --out ${CF_OUTPUT_FILE} --corr-type cf

#Export
${PICCA_PATH}/bin/picca_export.py --data ${CF_OUTPUT_FILE} --out ${CF_OUTPUT_EXP_FILE}
date

#Generate the aux files
/global/homes/j/jfarr/Projects/LyaCoLoRe/picca_analysis/make_picca_aux_files.py --base-dir ${OUTPUT_PATH} --corr-type cf --quantity $QUANTITY --npixels ${NPIXELS} --quant-code ${CF_QC} --zmin ${CF_ZMIN} --zmax ${CF_ZMAX} --cf-filename $CF_OUTPUT_FILENAME --cf-exp-filename $CF_OUTPUT_EXP_FILENAME --nside $CF_NSIDEPICCA --rpmin ${CF_RPMIN} --rpmax ${CF_RPMAX} --rtmin ${CF_RTMIN} --rtmax ${CF_RTMAX} --np ${CF_NP} --nt ${CF_NT} --rmin-values ${RMINS} --rmax-values ${RMAXS} --afix-values ${AFIXS}

EOF
fi

if [ $XCF -eq 1 ]; then
XCF_RUN_FILE="${OUTPUT_PATH}/run_picca_xcf_${ANALYSIS_ID}_${PICCA_ID}.sh"
if [ ! -z ${XCFSHUFFLESEED} ]; then NNODES_XCF_TOT=$(( 2 * $NNODES ))
else NNODES_XCF_TOT=$NNODES; fi
echo "xcf run file "${XCF_RUN_FILE}
cat > $XCF_RUN_FILE <<EOF
#!/bin/bash -l

#SBATCH --partition ${QUEUE}
#SBATCH --nodes ${NNODES_XCF_TOT}
#SBATCH --time ${TIME}
#SBATCH --job-name picca-xcf-${QC}-${NPIXELS}
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

    #Run xcf without shuffling.
    NODE_OUTPUT_FILE="${NODE_OUTPUT_PATH}/xcf_\${START_INDEX}_\${STOP_INDEX}.fits.gz"
    command="srun -N 1 -n 1 -c ${NCORES} ${PICCA_PATH}/bin/picca_xcf.py --in-dir \"dummy\" --from-image \${NODE_FILES} --drq ${DRQ} --out \${NODE_OUTPUT_FILE} --rp-min ${XCF_RPMIN} --rp-max ${XCF_RPMAX} --rt-max ${XCF_RTMAX} --np ${XCF_NP} --nt ${XCF_NT} --no-project --no-remove-mean-lambda-obs --nside ${XCF_NSIDEPICCA} --nproc 64 --z-cut-min ${XCF_ZMIN} --z-cut-max ${XCF_ZMAX} --z-min-obj ${XCF_ZMINQSO} --z-max-obj ${XCF_ZMAXQSO}";
    echo \$command
    \$command >& ${OUTPUT_PATH}/logs/node-\${NODE}.log &

    #If desired, run xcf with shuffling.
    if [ ! -z ${XCFSHUFFLESEED} ]; then

        NODE_OUTPUT_FILE="${NODE_OUTPUT_PATH}/xcf_shuffle_\${START_INDEX}_\${STOP_INDEX}.fits.gz"
        command="srun -N 1 -n 1 -c ${NCORES} ${PICCA_PATH}/bin/picca_xcf.py --in-dir \"dummy\" --from-image \${NODE_FILES} --drq ${DRQ} --out \${NODE_OUTPUT_FILE} --rp-min ${XCF_RPMIN} --rp-max ${XCF_RPMAX} --rt-max ${XCF_RTMAX} --np ${XCF_NP} --nt ${XCF_NT} --no-project --no-remove-mean-lambda-obs --nside ${XCF_NSIDEPICCA} --nproc 64 --z-cut-min ${XCF_ZMIN} --z-cut-max ${XCF_ZMAX} --z-min-obj ${XCF_ZMINQSO} --z-max-obj ${XCF_ZMAXQSO} --shuffle-distrib-obj-seed ${XCFSHUFFLESEED}";
        echo \$command
        \$command >& ${OUTPUT_PATH}/logs/node-\${NODE}-shuffle.log &

    fi

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
/global/homes/j/jfarr/Projects/LyaCoLoRe/picca_analysis/combine_picca_cf_files.py --in-dir ${NODE_OUTPUT_PATH} --out ${XCF_OUTPUT_FILE} --corr-type xcf

#If needed, combine the shuffled output files together
if [ ! -z ${XCFSHUFFLESEED} ]; then
/global/homes/j/jfarr/Projects/LyaCoLoRe/picca_analysis/combine_picca_cf_files.py --in-dir ${NODE_OUTPUT_PATH} --out ${XCF_SHUFFLE_OUTPUT_FILE} --corr-type xcf_shuffle
fi

#Export
if [ ! -z ${XCFSHUFFLESEED} ]; then
${PICCA_PATH}/bin/picca_export.py --data ${XCF_OUTPUT_FILE} --out ${XCF_OUTPUT_EXP_FILE} --remove-shuffled-correlation ${XCF_SHUFFLE_OUTPUT_FILE}
else
${PICCA_PATH}/bin/picca_export.py --data ${XCF_OUTPUT_FILE} --out ${XCF_OUTPUT_EXP_FILE}
fi
date

#Generate the aux files
/global/homes/j/jfarr/Projects/LyaCoLoRe/picca_analysis/make_picca_aux_files.py --base-dir ${OUTPUT_PATH} --corr-type xcf --quantity $QUANTITY --npixels ${NPIXELS} --quant-code ${XCF_QC} --zmin ${XCF_ZMIN} --zmax ${XCF_ZMAX} --cf-filename $XCF_OUTPUT_FILENAME --cf-exp-filename $XCF_OUTPUT_EXP_FILENAME --nside $XCF_NSIDEPICCA --rpmin ${XCF_RPMIN} --rpmax ${XCF_RPMAX} --rtmin ${XCF_RTMIN} --rtmax ${XCF_RTMAX} --np ${XCF_NP} --nt ${XCF_NT} --rmin-values ${RMINS} --rmax-values ${RMAXS} --afix-values ${AFIXS}

EOF
fi

# send the jobs to the queue
if [ $CF -eq 1 ]; then
sbatch $CF_RUN_FILE
fi
if [ $XCF -eq 1 ]; then
sbatch $XCF_RUN_FILE
fi
