# Umbrella_Sampling

This repository contains the data required and results of umbrella sampling of NaCl system obtained from MBAR analysis.

`data`: This folder contains all of the input files for MBAR, including one `center.dat` file and 47 `pullx-umbrella.xvg` files.  

`center.dat`: This file contains the data of the force constant used for each simulation.  

`pullx-umbrella.xvg`: There are 47 .xvg files which contain the data (the value of reaction coordinate as a functino of time) extracted from 47 independent simulations along the reaction coordinate.  

`NaCl_analysis_MBAR_bootstrapping.py`: This is the python code for MBAR analysis. PMF data can be obtained by simply running the file and all the figures will be save automatically, including `Bootstrapping_MBAR_var.png` and `PMF_MBAR.png`. The section of bootstrapping starts from line 168 and can be turned on or off in line 22. (The number of bootstrap samples can be defined in lie 23.)

`PMF_MBAR.png`: This figure was generated from `NaCl_analysis_MBAR_bootstrapping.py` with 200 bins used for plotting.  

`Bootstrapping_MBAR_var.png`: This figure is the result of bootstrapping of MBAR obtained from `NaCl_analysis_MBAR_bootstrapping.py`, with 200 bootstrap samples and error bars as the variances of samples. (However, the error bars are not obvious, since the variances are really small (around 0.02). 

`hist.xvg` and `profile.xvg`: These files were generated from the execution of g_wham command (`gmx_mpi wham -it ../../analysis/tpr_list.dat -if ../../analysis/pullf_list.dat -
b 500 -e 50000 -o ../../analysis/profile.xvg -hist ../../analysis/histo.xvg -unit kT`). The first 500 ps of simulation was truncated and 200 bins (defulat number) were used. In the `profile.xvg`, the unit of PMF is kT. (T=300K) 

`PMF_WHAM.jpg`: This figure was generated using the data of `profile.xvg`.

`Histogram.jpg`: This figure was generated using the data of `hist.xvg`.

`PMF_compare.jpg`: This figure shows a comparison between PMF curves obtained from PMF and MBAR.



