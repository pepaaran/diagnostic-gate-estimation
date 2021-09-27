# Code for diagnosis of interference regression adjustments
This repository contains code for producing the simulations and plots in my master thesis: Diagnostic Tool for Interference Features in GATE Estimation.

The code available at https://github.com/ajchin/regression-adjustments was used as a starting point for our implementations and simulations. 
The repository contains code used to produce the results shown in A. Chin, Regression adjustments for estimating the global treatment effect in experiments with interference, Journal of Causal Inference, May 2019 (https://www.degruyter.com/view/j/jci.ahead-of-print/jci-2018-0026/jci-2018-0026.xml)

## Data
The directory ``data`` serves as a placeholder folder for the data.

- ``facebook100`` is a placeholder for the files from the Facebook100 dataset, which is publicly available from the Internet Archive at https://archive.org/details/oxford-2005-facebook-matrix and other public repositories. 

## Reproducing results and plots
The directory ``src`` contains all the code used for the thesis report

- ``setup.R`` installs and loads all the packages used in the simulations
1. The  ``functions`` folder contains helpers classified as follows:
- ``covariate_functions.R`` contains functions for generating some possible interference features. 
- ``response_functions.R`` generates the outcomes for different models.
- ``data_generators.R`` builds data frames containing realized samples of the data.
- ``existing_estimators.R`` contains baseline estimators used in Chin (2019).
- ``proposed_estimators.R`` constructs the adjustment estimators from Chin (2019).
- ``variance_estimators.R`` computes the estimated variance for the proposed estimators.
- ``precompute_matrices.R`` contains an intermediate function that estimates the matrices involved in the variance estimation.
- ``permutation_tests.R`` contains our implementation for the proposed interference tests
2. The ``plots`` directory contains the code used to generate the plots in the thesis report.
3. The ``simulations_final.R`` directory contains all the scripts and config files used in simulations (run on the ETH Euler cluster)

 The ``results_final`` directory is a placeholder for the simulation outputs, stored as .csv files.
