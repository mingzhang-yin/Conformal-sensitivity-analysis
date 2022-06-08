# Conformal Sensitivity Analysis 

Implemention of conformal sensitivity analysis (CSA) described in [Conformal Sensitivity Analysis for Individual Treatment Effects](https://arxiv.org/pdf/2112.03493.pdf) by Mingzhang Yin, Claudia Shi, Yixin Wang and David Blei

## This repository

* Provides code to reproduce the simulation results in the paper
* Implements sensitivity analysis for the ITE on synthetic data and an observational study 
* Requires the [cfcausal package](https://github.com/lihualei71/cfcausal), which can be installed according to the author's instructions
```
library("devtools")
library("cfcausal")
devtools::load_all(".")
```
## Directories

* The Conformal Sensitivity Analysis (CSA) algorithm is implemented by the R scripts in the folder `./R`
* The Conformalized Sharp Sensitivity Analysis (CSSA) algorithm is in `./exp-cf/synthetic_exp-sharp.R`
* The folders `./exp-cf`, `./exp-ite`, `./exp-fish` are the example folders that contain the executable scripts for the simulations in the paper

## Code scripts, results and data

* Each example folder has a directory `./bash` that contain the bash scripts. The bash scripts are used to submit jobs to a cluster, call the corresponding R scripts, and  produce the results 
```
cd ./bash
bash run_XXX.sh
```
* The python script in each example folder can be used to generate the plots in the paper
* Each example folder has a directory `./out` that contains current numerical results 
* The fish consumption and blood mercury data is from R package [bootsens](https://github.com/qingyuanzhao/bootsens). The data is also contained in `./exp-fish/data`







