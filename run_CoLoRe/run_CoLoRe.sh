# specify grid size
NGRID=4096

# specify number of nodes to use
NODES=32

# specify seed
SEED=$1

# specify queue
QUEUE=$2

# smoothing radius
R_SMOOTH=2.0

# gaussian skewers?
GAUSSIAN="true"

# full path to CoLoRe_revamp executable (parallel version)
COLORE_PATH="/global/homes/j/jfarr/Programs/CoLoRe_PAR_GAUSS/"

# full path to folder where output will be written
OUTPUT_PATH=$3
echo "output will written to "$OUTPUT_PATH
mkdir $OUTPUT_PATH

# if it does not exist already, make a run_files directory
if [ ! -d $OUTPUT_PATH/run_files ] ; then
    mkdir -p $OUTPUT_PATH/run_files
fi

# we will create this parameter file
PARAM_FILE="${LYACOLORE_PATH}/run_CoLoRe/run_files/param_seed${SEED}_${NGRID}.cfg"
echo "parameter file "$PARAM_FILE

cat > $PARAM_FILE <<EOF
global = {
  prefix_out = "${OUTPUT_PATH}/out"
  output_format = "FITS"
  output_density = false
  pk_filename = "${LYACOLORE_PATH}/run_CoLoRe/input_files/PlanckDR12_kmax_matterpower_z0.dat"
  z_min = 1.6
  z_max = 3.79
  seed = ${SEED}
  write_pred = false
  pred_dz = 0.1
}
field_par = {
  r_smooth = $R_SMOOTH
  smooth_potential = true
  n_grid = $NGRID
  dens_type = 0
  lpt_buffer_fraction = 0.6
  lpt_interp_type = 1
  output_lpt = 0
}
cosmo_par = {
  omega_M = 0.3147
  omega_L = 0.6853
  omega_B = 0.04904
  h = 0.6731
  w = -1.0
  ns = 0.9655
  sigma_8 = 0.830
}
srcs1 = {
  nz_filename = "${LYACOLORE_PATH}/run_CoLoRe/input_files/Nz_qso_130618_2_colore1_hZs.txt"
  bias_filename = "${LYACOLORE_PATH}/run_CoLoRe/input_files/Bz_qso_G18.txt"
  include_shear = false
  store_skewers = true
  gaussian_skewers = ${GAUSSIAN}
}

EOF

# we copy the param file to the output directory
cp $PARAM_FILE $OUTPUT_PATH

# we will create this script
RUN_FILE="${LYACOLORE_PATH}/run_CoLoRe/run_files/run_colore_seed${SEED}_${NGRID}.sh"
echo "run file "$RUN_FILE

cat > $RUN_FILE <<EOF
#!/bin/bash -l

#SBATCH --partition $QUEUE
#SBATCH --nodes $NODES
#SBATCH --time 00:30:00
#SBATCH --job-name CoLoRe_Lya_${NGRID}_${NODES}_highZ
#SBATCH --error "${LYACOLORE_PATH}/run_CoLoRe/run_files/CoLoRe_Lya_${NGRID}_${NODES}-%j.err"
#SBATCH --output "${LYACOLORE_PATH}/run_CoLoRe/run_files/CoLoRe_Lya_${NGRID}_${NODES}-%j.out"
#SBATCH -C haswell
#SBATCH -A desi

umask 0002
export OMP_NUM_THREADS=64

date
srun -n $NODES -c 64 ${COLORE_PATH}/CoLoRe $PARAM_FILE  
date

EOF

sbatch $RUN_FILE

# we copy the param file to the output directory
cp $RUN_FILE $OUTPUT_PATH
