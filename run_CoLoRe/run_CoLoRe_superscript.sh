# specify the number of boxes to run
NBOXES=1

# specify which seed to start at (subsequent seeds chosen sequentially)
START_SEED=1003
SEEDS=$(seq $START_SEED $((START_SEED + NBOXES - 1)))

#specify which queue to use
#<5 boxes can use debug, otherwise regular
QUEUE='debug'

for SEED in $SEEDS; do

# specify where to output the box to
OUTDIR="/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/CoLoRe_seed${SEED}_4096/"

# run colore
source run_CoLoRe.sh ${SEED} ${QUEUE} ${OUTDIR};

done
