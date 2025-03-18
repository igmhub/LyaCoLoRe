COLORE_BOX CONTAINS FILES TO REPRODUCE COLORE_QL BOXES

- The CoLoRe QL boxes used the linear biasing model with a threshold of CoLoRe (bias model 3 in the Makefile)

INPUT_FILES FOLDER CONTAINS:
- bias_colore-ql.txt, Nz.txt and threshold_colore-ql.txt  --> 1st column all of these .txt is the redshift. 2nd column is the input bias, input dN/dzdOmega and input threshold respectively.
    
- smoothed_mu0_NL_matter_powerspectrum_DR12.dat is the input power spectrum at z=0 used to generate the CoLoRe-QL boxes. It has the smoothing of the BAO peak added for mu=0 and uses DR12 cosmology.
    
    
NOTEBOOKS FOLDER CONTAINS:
    
- run_colore-ql_box.ipynb will send the job to run the colore-ql boxes. Note that the random seed will be 2024 + number of the box. Hence box 0 will have random seed 2024. 




LYACOLORE_SKEWERS

- lyacolore_run_box.sh --> script to run the lyacolore skewers from a colore box. Run on terminal doing bash /lyacolore_run_box.sh
- config_v9.0. --> config file to point to on the lyacolore_run_box.sh