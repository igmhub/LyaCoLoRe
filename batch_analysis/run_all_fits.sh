BASEDIR="/project/projectdirs/desi/users/jfarr/job_submit_tests/test_1/"
cd $BASEDIR
CHI2_FILES=`ls analysis/correlation_functions/v9.0.0/fits/lya_dla_cross/chi2*.ini`
echo $CHI2_FILES
for CHI2_FILE in $CHI2_FILES; do
/global/homes/j/jfarr/Programs/picca/bin/picca_fitter2.py ${CHI2_FILE}
done

cd -

