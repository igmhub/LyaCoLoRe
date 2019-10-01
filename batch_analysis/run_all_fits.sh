BASEDIR="/project/projectdirs/desi/users/jfarr/LyaCoLoRe_paper/"
cd $BASEDIR

# Set up the versions to cycle through.
RS=`seq 0 9`
VERS=""
for R in $RS; do
	VERS="$VERS v9.0.$R"
done
VERS="$VERS stack"

# Set up the measurements to cycle through.
#CTS="lya_auto dla_auto lya_qso_cross lya_dla_cross lya_auto__lya_qso_cross lya_auto__lya_dla_cross qso_auto lya_aa_auto"
CTS="dla_auto qso_auto lya_auto__lya_dla_cross__dla_auto lya_auto__lya_qso_cross__qso_auto"
#CTS="lya_auto__lya_dla_cross__dla_auto"

# Cycle through everything.
for VER in $VERS; do
	echo $VER
	for CT in $CTS; do
		echo $CT
		CHI2_FILES=`ls analysis/correlation_functions/$VER/fits/$CT/chi2*.ini`
		for CHI2_FILE in $CHI2_FILES; do
			/global/homes/j/jfarr/Programs/picca/bin/picca_fitter2.py ${CHI2_FILE}
	
		done
	done
done

cd -

