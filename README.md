# Code for the paper "Doubly robust, machine learning effect estimation in real-world clinical sciences: A practical evaluation of performance in molecular epidemiology cohort settings"

Code to reproduce results in the paper [Doubly robust, machine learning effect estimation in real-world clinical sciences: A practical evaluation of performance in molecular epidemiology cohort settings](https://arxiv.org/abs/2105.13148) by Xiang Meng, Jonathan Huang

Xiang Meng 2021.06.11

## 1. Code Structure
* [2020309-Plasmode_Sim_CV-TMLE.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/2020309-Plasmode_Sim_CV-TMLE.R): Master code to first generate dataset and then produce estimation. It has dependencies on all other codes. 
* [20200803-Sims-Function.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/20200803-Sims-Function.R): Code to generate simulation dataset and plasmode dataset
* [20200226-PlasmodeCont_Revised.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/20200226-PlasmodeCont_Revised.R): Function to generate Plasmode dataset
* [20200904-run-sim-code.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/20200904-run-sim-code.R): Code to parallelize estimation. 
* [20200720-Algos-code.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/20200720-Algos-code.R): Function that contains all estimators
* [20200705-DCDR-Functions.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/20200705-DCDR-Functions.R): Function for Double Cross-Fit which is modified from Breskin and Zivich (2020).
* [20200816-Result-Summary.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/20200816-Result-Summary.R): A script containing all summary functions.
* [3.1 Single-Summary.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/3.1%20Single-Summary.R): Code that produces summary for a single result file (need to manually change path)
* [3.Summary.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/3.Summary.R): Code that produces a list of summaries and plots for the paper.
* [OutputCoeff.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/OutputCoeff.R): Code to produce the actual simulation coefficients as shown in the appendix.

## 2. Results in Paper and Corresponding Code
As briefly mentioned above, parameters in [2020309-Plasmode_Sim_CV-TMLE.R](https://github.com/mengeks/drml-plasmode/blob/master/Code/2020309-Plasmode_Sim_CV-TMLE.R) controls all the simulations. To run for different simulation scenario, one should just change parameters here. eg. "idx.handpick" is the indices of covariates picked; "interact" means whether we add first order interactions between covariates; "interact.w.exp" means we add first order interactions between covariates and the exposure variable. 
