# LyaCoLoRe
Code development to use CoLoRe simulations for generating simulated Lyman alpha forest spectra

## Install
You only need to include LyaCoLoRe/py in your PYTHONPATH, with something like: 
export PYTHONPATH=$PYTHONPATH:$HOME/Programs/igmhub/LyaCoLoRe/py

## Examples
You can find some examples under example_scripts/. Specifically, the examples available are:

 - plot_delta_picca.ipynb
      To quickly open and view files in the picca input format.

 - plot_raw_colore.ipynb
      To quickly open and view files in the CoLoRe output format.

 - mock_desi.ipynb (requires DESI specific packages)
      To produce mock DESI skewers.

 - read_skewers_lya.ipynb
      To look in detail at files in the CoLoRe output format.

In addition, issue_dec_colore.ipynb and issue_mean_delta_colore.ipynb describe two ongoing issues with CoLoRe.
