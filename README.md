# cvtmle_plasmode
Diagnostics of doubly-robust estimators in real world (plasmode-simulated) finite samples

Current best practices for quantitative causal inference recommend machine learning (non-parametric/non-smooth classifers) of nuisance parameters be used with care: namely, using doubly robust estimators (treatment and outcome models; Naimi, Mishler, Kennedy 2020) and cross-fitting (parameter and standard error estimation on a separate set than the one used for nuisance function estimation; Zivich & Breskin 2020). Moreover, despite reasonably faster convergence rates for these doubly robust estimators, small sample sizes may be insufficient to guarantee appropriate coverage properties regardless of the estimators implemented (Benkeser 2017, Rotnitzky 2019). 

Most importantly, existing simulation studies use fairly simple confounding structures (Naimi 2020) or large sample sizes (Bahamyirou 2018, Zivich 2020) which may demonstrate overly optimistic performance relative to the real world use case of complex confounding, practical positivity violations, and small samples. Here we demonstrate performance in closer to real-world conditions, notably by applying estimators to a structure-preserving "plasmode" dataset (where covariate data are retained but treatment and outcome values reassigned to force a known parameter estimate across bootstrapped data sets) drawn from an existing longitudinal cohort study (N = 1178; 331 covariates). One past effort by Bahamyirou, et al. (2018) focusing on positivity violations in propensity score estimation utilized a similar approach (Bootstrap Diagnostic Test) but on large samples an fewer, simple covariates.

Overall objectives:
1. Reproduce earlier findings showing that optimal bias and coverage is provided by cross-fit DR (TMLE & AIPW) estimators when the data sets are large (N > 2000) and the confounding sets are simulated simply (with or without misspecfication)
2. Show that this does not hold when applied to "real world" data which contain more complex covariation and likelihood for positivity violations using "plasmode" (structure-preserving) simulation based on an available cohort dataset (N <= 1178; 16 or 331 covariates)
3. Provide some guidance on diagnotics (using plasmode) and best practices in these real-world sets

What has been done:
1. Function to sample plasmode data sets was created based on existing (non-functional) "plasmode" package (0200226-PlasmodeCont_Revised.R) 
2. Simulations demonstrating good performance of TMLE / AIPW in an "ideal setting" (e.g. 5 covariates including non-linearities, N = 3000, 10k bootstrap samples); code draws simulations, fits AIPW, CV-TMLE, and IPW+GLM models, estimates contrasts, computes bias and coverage (20200304-Basic_Sim_CV-TMLE.R)
3. Some simulations for the "real world" setting (16 or 331 covariates; N = 1178) varying the number of bootstrapped sets (200-1000) and fitting AIPW, CV-TMLE, and IPW+GLM using smooth (glm, elastic net, polynomial splines) and non-smooth (random forest, xgboost, single-layer neural net) estimators of nuisance parameters, computing bias and coverage (20200309-Plasmode_Sim_CV-TMLE.R)  

Preliminary findings:
We reproduce the finding that cross-fit, doubly robust estimators do have reduced bias and improved coverage to simple regression in simple confounding cases and larger datasets. In plasmode simulated data with slightly more covariates (N = 16) cross-fit AIPW and TMLE work well with nuisance parameters estimated with smooth functions (glm, elastic net, polynomial splines) but, more importantly, are biased with poor coverage if non-parametric classifiers are used (e.g. random forest, xgboost). This was true even when following best-practices for a moderate sized sample (N = 1178). A simple IPW+GLM performed just as well. In the case of high-dimensional covariates, where some variable selection / penalization is necessary (i.e. machine learning) non-smooth nuisance estimation was able to recover the target contrast, but coverage remained suboptimal.  

To-Do (as of 9 June 2020):
1. Verify code for plasmode simulation, "ideal setting," and "real world" setting performance
2. Optimize the CV-AIPW in some way (otherwise ensemble learner must be fit 3 times per run)
3. Complete simulations for 200 and 500 samples (perhaps do not need to vary the number of bootstraps used to estimate -- not much variation in this regard)
4. Complete interpretation and write-up
