################################################################################
BASEDIR="/project/projectdirs/desi/users/jfarr/LyaCoLoRe_paper/"
LYACOLORE_PATH="/global/homes/j/jfarr/Projects/LyaCoLore/"
RANDDIR="${BASEDIR}/data/additional_data/"
V_CODE_MAJOR=9
V_CODE_MINOR=0
V_REALISATIONS=`echo {0..9}`
NPROC=64
NSIDE=16
FLAGS="--make-zcats --make-random-zcats"
#--DLAs-in-transmission-rest-range --single-DLA-per-skw --add-Lyb --add-metals

#Processing quantities.
DS=0.5
DLA_Z_BUFFER=0.05
DS_RAND=0.1
#Set the seeds for the downsampling using the version number.
MIN_CAT_Z=1.7

# Wavelength grid for output.
lObs_min=3600.
lObs_max=5500.
lRF_min=1040.
lRF_max=1200.
dll=0.0003

################################################################################
## For each realisation, run script to make the deltas, once using only lya
## absorption, once using all absorbers.
for r in $V_REALISATIONS; do
INDIRNAME="$BASEDIR/data/LyaCoLoRe_output/v${V_CODE_MAJOR}.${V_CODE_MINOR}.${r}/"
OUTDIRNAME="$BASEDIR/data/picca_input/v${V_CODE_MAJOR}.${V_CODE_MINOR}.${r}/"
RUN_FILE="$OUTDIRNAME/run_make_deltas_v${V_CODE_MAJOR}.${V_CODE_MINOR}.${r}.sh"
if [ ! -d $OUTDIRNAME ] ; then
    mkdir -p $OUTDIRNAME
fi

################################################################################
RUN_FILE="$OUTDIRNAME/run_make_deltas_v${V_CODE_MAJOR}.${V_CODE_MINOR}.${r}.sh"
echo "Run file for making picca deltas (lya only) will be written to "$RUN_FILE
cat > $RUN_FILE <<EOF
#!/bin/bash -l

#SBATCH --partition debug
#SBATCH --nodes 1
#SBATCH --time 00:30:00
#SBATCH --job-name run_make_deltas_v${V_CODE_MAJOR}.${V_CODE_MINOR}
#SBATCH --error "${OUTDIRNAME}/run-make-deltas-%j.err"
#SBATCH --output "${OUTDIRNAME}/run-make-deltas-%j.out"
#SBATCH -C haswell
#SBATCH -A desi

umask 0002
export OMP_NUM_THREADS=64

python $LYACOLORE_PATH/scripts/make_picca_deltas_from_transmission.py --in-dir $INDIRNAME --out-dir $OUTDIRNAME --randoms-dir $RANDDIR --nproc $NPROC --downsampling $DS --downsampling-seed $r --randoms-downsamping $DS_RAND --randoms-downsampling-seed $(( $r + 1000 )) --min-cat-z $MIN_CAT_Z --lambda-obs-min $lObs_min --lambda-obs-max $lObs_max --lamdba-RF-min $lRF_min --lambda-RF-max $lRF_max $FLAGS

wait
date

EOF

################################################################################
## Send the job to the queue.
sbatch $RUN_FILE

################################################################################

RUN_FILE="$OUTDIRNAME/run_make_deltas_v${V_CODE_MAJOR}.${V_CODE_MINOR}.${r}_all_absorbers.sh"
echo "Run file for making picca deltas (all absorbers) will be written to "$RUN_FILE
cat > $RUN_FILE <<EOF
#!/bin/bash -l

#SBATCH --partition debug
#SBATCH --nodes 1
#SBATCH --time 00:30:00
#SBATCH --job-name run_make_deltas_v${V_CODE_MAJOR}.${V_CODE_MINOR}
#SBATCH --error "${OUTDIRNAME}/run-make-deltas-%j.err"
#SBATCH --output "${OUTDIRNAME}/run-make-deltas-%j.out"
#SBATCH -C haswell
#SBATCH -A desi

umask 0002
export OMP_NUM_THREADS=64

python $LYACOLORE_PATH/scripts/make_picca_deltas_from_transmission.py --in-dir $INDIRNAME --out-dir $OUTDIRNAME --randoms-dir $RANDDIR --nproc $NPROC --downsampling $DS --downsampling-seed $r --randoms-downsamping $DS_RAND --randoms-downsampling-seed $(( $r + 1000 )) --min-cat-z $MIN_CAT_Z --lambda-obs-min $lObs_min --lambda-obs-max $lObs_max --lamdba-RF-min $lRF_min --lambda-RF-max $lRF_max $FLAGS --add-Lyb --add-metals

wait
date

EOF

done

################################################################################
## Send the job to the queue.
sbatch $RUN_FILE

################################################################################ 
