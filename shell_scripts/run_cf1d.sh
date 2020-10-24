#skewers or revamp
NPROCESSES=64

#branch
BRANCH='revamp'

#option
OPTION='gaussian'

# we will create this script
RUN_FILE="/global/homes/j/jfarr/Projects/LyaCoLoRe/run_files/run_cf1d.sh"
echo "run file "$RUN_FILE

cat > $RUN_FILE <<EOF
#!/bin/bash -l

#SBATCH --partition regular
#SBATCH --nodes 1
#SBATCH --time 02:00:00
#SBATCH --job-name cf1d
#SBATCH --error "/global/homes/j/jfarr/Projects/LyaCoLoRe/run_files/run-cf1d-%j.err"
#SBATCH --output "/global/homes/j/jfarr/Projects/LyaCoLoRe/run_files/run-cf1d-%j.out"
#SBATCH -C haswell
#SBATCH -A desi

umask 0002
export OMP_NUM_THREADS=64

date

/global/homes/j/jfarr/Projects/LyaCoLoRe/example_scripts/cf1d.py ${NPROCESSES}

date

EOF

sbatch $RUN_FILE
