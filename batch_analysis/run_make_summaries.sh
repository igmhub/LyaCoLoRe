BASEDIR="/project/projectdirs/desi/users/jfarr/LyaCoLoRe_paper/"
LYACOLORE_PATH="/global/homes/j/jfarr/Projects/LyaCoLore/"
V_CODE_MAJOR=9
V_CODE_MINOR=0
V_REALISATIONS=`echo {0..9}`
NPROC=32
NSIDE=16
FLAGS="--compressed-input --compress --transmission-only"

## For each realisation, run make_summaries.py
for r in $V_REALISATIONS; do
DIRNAME="$BASEDIR/data/v${V_CODE_MAJOR}.${V_CODE_MINOR}.${r}/LyaCoLoRe_output/
python $LYACOLORE_PATH/scripts/make_summaries.py --base-dir $DIRNAME --nproc $NPROC --nside $NSIDE $FLAGS
done
