# DHMC for Stochastic-Daily-Rainfall-Simulation
This repository presents a Python code of the Decadal and Hierarchical Markov Chain (DHMC) model for stochastic daily rainfall simulation. Technical details of the model are available in the following publications. 

Chowdhury, A.F.M.K., Lockart, N., Willgoose, G., Kuczera, G., Kiem, A. S., and Parana Manage, N. (2017). Development and evaluation of a stochastic daily rainfall model with long-term variability. Hydrology and Earth System Sciences, 21(12), 6541-6558. DOI: [10.5194/hess-21-6541-2017](https://www.hydrol-earth-syst-sci.net/21/6541/2017/hess-21-6541-2017.html).

Chowdhury, A.F.M.K. (2017). Development and evaluation of stochastic rainfall models for urban drought security assessment, PhD Thesis, The University of Newcastle, Australia. [http://hdl.handle.net/1959.13/1342370](http://hdl.handle.net/1959.13/1342370)

**If you use the code, please make sure to cite the above publications.**



# Contents

**data_site_$$.txt:** sample data of observed daily rainfall of two test sites.

**DHMC_Calibration_Simulation.py:** This script calibrates the parameters of DHMC and generates multiple (say 1000) replicates of stochastic daily rainfall series. Calibrated parameters and simulated rainfall can be saved in text files. 

**Calc_ZScore_from_SimRain.py:** As a sample test of model performance, this script calculates Z Scores for rainfall depth at daily, monthly and multiyear resolutions and saves in textfiles.
