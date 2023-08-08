# EP-ABC

This repository contains the code for running `Expectation-Propagation-ABC` for a Neutral evolutionary model. 

We present a new model-based simulation approach to study evolution and demography from whole-genome data. It incorporates both coalescent and forward simulations and parallel Expectation Propagation ABC (EP-ABC). EP-ABC makes use of recent approaches in machine-learning and a “divide to conquer” approach to enable efficient distributed computation in genomic analysis. It uses a a population size history model where we  allow effective population size fluctuations through time.



There is a `Template folder` that contains all the files needed to run the `EP-ABC` and an a `CreateInfiles_onepop.sh`. The user needs to change the parameters on this file and run it locally to generate the infiles.

It runs in the `PBS Pro` job scheduler for HPC. Changes need to be done to the `bash` files to converting to them `Slurm`.