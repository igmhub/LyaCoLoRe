BASEDIR="/project/projectdirs/desi/users/jfarr/LyaCoLoRe_paper/"
LYACOLORE_PATH="/global/homes/j/jfarr/Projects/LyaCoLore/"
RANDDIR="${BASEDIR}/additional_data/"
V_CODE_MAJOR=9
V_CODE_MINOR=0
V_REALISATIONS=`echo {0..9}`
NPROC=32
NSIDE=16
FLAGS="--compressed-input --compress --add-Lyb --add-metals --downsample-randoms"

#Processing quantities.
DS=0.5
DS_SEED=42
DLA_Z_BUFFER=0.05

# Wavelength grid for output.
lObs_min=3600.
lObs_max=5500.
lRF_min=1040.
lRF_max=1200.
dll=0.0003


## For each realisation, run make_summaries.py
for r in $V_REALISATIONS; do
INDIRNAME="$BASEDIR/data/LyaCoLoRe_output/v${V_CODE_MAJOR}.${V_CODE_MINOR}.${r}/"
OUTDIRNAME="$BASEDIR/data/picca_input/v${V_CODE_MAJOR}.${V_CODE_MINOR}.${r}/"
RUN_FILE="$OUTDIRNAME/run_make_deltas_v${V_CODE_MAJOR}.${V_CODE_MINOR}.${r}.sh"

echo "Run file for making picca deltas will be written to "$RUN_FILE
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

python $LYACOLORE_PATH/scripts/make_picca_deltas_from_transmission.py --in-dir $INDIRNAME --out-dir $OUTDIRNAME --randoms-dir $RANDDIR --nproc $NPROC --downsampling $DS --downsampling-seed $DS_SEED --DLA-z-buffer $DLA_Z_BUFFER --lambda-obs-min $lObs_min --lambda-obs-max $lObs_max --lamdba-RF-min $lRF_min --lambda-RF-max $lRF_max $FLAGS

wait
date

EOF

done

