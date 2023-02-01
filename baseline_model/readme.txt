####### Predicting Food Crises Code and Data Files ########
# Written by: Bo Pieter Johannes Andree 
# Last edited: 08.10.2020


CONTENTS OF THIS FILE
---------------------

 * Introduction
 * Requirements
 * Installation
 * Configuration
 * Execution
 * Troubleshooting
 * Maintainers
 * sessionInfio

INTRODUCTION
------------

This package contains code and data for a statistical foresting approach to predict the outbreak of food crises. 
Different use cases are explored related to possible alternative targeting policies and the levels at which finance is typically unlocked. 

The package consists of the following files:

 * predicting_food_crises_data.csv
     - data set used to produce results of the paper. 
 * metadata.xlsx
     - a description of the data set and sources. 
 * predicting_food_crises.R
     - the main R code used to produce to process data, train and cross-validate models for 1 to 12 month ahead forecasts and generate useful plots and results.
 * predicting_food_crises_dependencies.R
     - dependencies required to run the main R code contained in predicting_food_crises.R
 * predicting_food_crises_balanced_learners.R
     - popular classifiers with cross-validated probability balancing as described in the paper. 
     - the code can be run using caret by supplying the object names as the "method=" argument in caret's train function and can be used for other prediction problems.
 * Predicting-Food-Crises.pdf
     - The Policy Research Working Paper (WPS 9412). It can also be found at the following URI: http://hdl.handle.net/10986/34510
     - The citation of this work is Andree, Bo Pieter Johannes; Chamorro, Andres; Kraay, Aart; Spencer, Phoebe; Wang, Dieter. 2020. Predicting Food Crises. Policy Research Working Paper; No. 9412. World Bank, Washington, DC. 


REQUIREMENTS
------------

Reading files:
 * The data set can be read and viewed in any text editor or excel. 
 * The meta data can be read with any program that supports excel files.
 * The .R files can be read with any text editor, the code was written in Sublime Text 3.

Execution:
 * The code was developed and last ran in Microsoft Open R 3.5.1, on Ubuntu 16.04.5 LTS and has not been tested on other OS. It should run in R 3.5.1 but the code benefits from multithreaded BLAS/LAPACK and contains a call to automatically sets MKL threads. This may throw an error in R 3.5.1 but should not break the remainder of the code.
 * The results presented in the paper have been generated on a virtual machine with 64 CPUs and 256GB RAM. Producing the full set of results in the paper consumed around 12,000 core hours and was run on a D64s_v3 VM with 64CPUs and 256 GiB RAM. Some simplifications have been made to make the final code more usable, comments are left in the main R file. See also CONFIGURATION section. 

Viewing plots:
 * R Studio server is recommended, see INSTALLATION notes.


INSTALLATION
------------

The user will need to follow standard installation instructions for R. 
 * To avoid unexpected issues, it is recommended to run this code on a similar R installation and OS, i.e. Microsoft Open R 3.5.1. on Ubuntu 16.04.5 and r-studio-server 1.2.5001.

Install the required R packages (lines 5 - 34 in predicting_food_crises_dependencies.R). 
 * Note that many R packages require the user to install dependencies on ubuntu OS itself.
 * User will need to install packages manually, since currently, there is no good way to automatize this. This is due to the large number of (in)direct dependencies in and outside R.
   - At the end of this readme file, a print out of sessionInfo() is provided such that versions of all packages can be viewed.
 * Note that the main R code (predicting_food_crises.R) sources the dependencies, the balanced learners, and reads the data. 
   - The user needs to specify the folder that contains these files in line 8. The default value is:
     "/home/predicting _food_crisis_package/" 
     which assumes this package is unzipped in the home folder of ubuntu.

The code can be run in a terminal, in which case the data plots will not be visible to the user.
 * One solution is to run the code on R Studio server. When set up correctly, one can access the RStudio IDE from anywhere via a web browser and use plot functionality. The code was developed on r-studio-server 1.2.5001. This can be isntalled by following standard installation procedures.


CONFIGURATION
------------
There are a number of choices that the user can make to control the behavior of the main program:

 * Lines 15-26 are options to control the definition of the dependent variable and the treatment of independent variables. 
   The default settings runs a model on all countries, using ipc 3 and above as positive class, uses only exogenous covariates as predictors, 
   adds synthetic cases to the training data, calculates additional features, and restricts linear correlation to .75. These are the settings that correspond to the paper.
 * Lines 31-32 control the type of learner used, default settings correspond to a simplified RF algorithm that delivers good results (nearly identical to the paper) but runs much faster.
   See also the comments in the code.
 * Lines 35-37 control an imputation strategy in case a missing value is encountered, settings should not matter when the supplied data is used.
 * Lines 40-43 control the cross-validation, note that repetitions have been reduced to make the runtime and RAM requirements more manageable.
 * Lines 46-55 control the compute environment.

Default settings:
 * Note that parallel processing works differently on ubuntu than on other OS, but generally it involves generating copies of dependencies or compute environments and so memory requirements can be extremely high even when the initial data set seems manageable. For this reason the following simplifications have been made to default seetings:  
   - The number of validation samples has been reduced from 50 to 10.
   - The tuning parameters of the default RF model have been fixed at recommended values. To run full tuning or use one of the alternative balanced classifiers, change MODEL_METHOD to one of the classifiers from predicting_food_crises_balanced_learners.R
   - When an alterantive model is used, the length of the tuning grid has been reduced to 5, the paper uses 10.
 * These settings produce similar results as those presented in the main paper, but the runtime and RAM requirements have been drastically reduced (depending of course on the number of CPUs available). 
   - The final code at (recommended) default settings was last run on a D32s_v3 VM with 32CPUs and 128 GiB RAM, reaching 100% CPU utilization and approx 60% RAM utilization, and took just below 2.5 hours to complete with approximately another hour for the additional validation results.
   - By default, the code runs on the entire data set that is provided. Note that the paper only trains and cross-validates on data up to February 2019. With the current settings, it is thus straightforward to update the data set and make real forecasts.


EXECUTION
------------

Running code:
 * After installation, simply unpack the folder, point the code (line 8) to the correct folder and run predicting_food_crises.R.
 * The code is currently not set up to write results to disk. As always, complex R objects can be saved for re-use using saveRDS() and text can be written using write.csv().


TROUBLESHOOTING
------------

Dependencies:
 * Make sure all OS dependencies are installed such that all libraries can be installed. Then make sure that all R libraries are installed and that also their dependencies are installed.
 * See the sessionInfo() readout at the end of this file. 
 * Make sure the predicting_food_crises_dependencies.R and predicting_food_crises_balanced_learners.R files are correctly sourced.

Unexpected crash with different compute settings:
 * If a different VM is used or if changes are made to the settings and the program crashes halfway, then keep an eye on the RAM usage. On ubuntu this can be monitored using 
   > htop
   If RAM usage is too high, reduce the number of cores used in lines 46-55.

NA values in validation metrics:
 * A common issue with caret is that validation metrics return as NA. This is likely result of a missing dependency in the slave environment, which may occur because different OS handle parallelization differently. See if the issue persists when setting MODEL_METHOD to another value, for example "multinom".


MAINTAINERS
------------
This package does not receive long term support. For technical questions on the paper and code, reach out to:

Bo P.J. Andree at bandree(at)worldbank.org

For data-related questions, reach out to

Andres Chamorro at achamorroelizond(at)worldbank.org

For other questions, reach out to
Nadia Piffaretti at npiffaretti(at)worldbank.org





sessionInfio
------------

> sessionInfo()
R version 3.5.1 (2018-07-02)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.5 LTS

Matrix products: default
BLAS: /opt/microsoft/ropen/3.5.1/lib64/R/lib/libRblas.so
LAPACK: /opt/microsoft/ropen/3.5.1/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods
[8] base

other attached packages:
 [1] DescTools_0.99.24    Ecdat_0.3-1          Ecfun_0.1-7
 [4] fpp_0.5              tseries_0.10-45      lmtest_0.9-36
 [7] expsmooth_2.3        fma_2.3              xts_0.11-0
[10] psych_1.8.4          janitor_1.1.1        miceadds_2.13-63
[13] CDM_6.4-23           mvtnorm_1.0-8        doFuture_0.6.0
[16] future_1.9.0         doParallel_1.0.11    doSNOW_1.0.16
[19] snow_0.4-2           iterators_1.0.10     caret_6.0-86
[22] ggplot2_3.0.0        foreach_1.4.4        gplots_3.0.1
[25] dplyr_0.8.5          plyr_1.8.4           imputeTS_2.7
[28] phylin_1.1.1         mice_3.3.0           lattice_0.20-35
[31] MASS_7.3-50          pROC_1.12.1          spdep_0.7-7
[34] spData_0.2.9.0       Matrix_1.2-14        sp_1.3-1
[37] zoo_1.8-3            KRLS_1.0-0           timeSeries_3042.102
[40] timeDate_3043.102    forecast_8.4         TTR_0.23-3
[43] RevoUtils_11.0.1     RevoUtilsMath_11.0.0

loaded via a namespace (and not attached):
 [1] backports_1.1.2       GPArotation_2014.11-1 lazyeval_0.2.1
 [4] splines_3.5.1         polycor_0.7-9         listenv_0.7.0
 [7] fda_2.4.8             digest_0.6.15         gdata_2.18.0
[10] magrittr_1.5          sirt_2.7-50           cluster_2.0.7-1
[13] sfsmisc_1.1-2         recipes_0.1.10        globals_0.12.1
[16] gower_0.1.2           gmodels_2.18.1        jpeg_0.1-8
[19] colorspace_1.3-2      mitools_2.3           pan_1.6
[22] crayon_1.3.4          lme4_1.1-17           survival_2.42-6
[25] glue_1.3.0            gtable_0.2.0          ipred_0.9-6
[28] mirt_1.28             quantmod_0.4-13       dcurver_0.9.1
[31] jomo_2.6-2            scales_0.5.0          stinepack_1.4
[34] Rcpp_1.0.4.6          foreign_0.8-71        stats4_3.5.1
[37] lava_1.6.2            survey_3.33-2         prodlim_2018.04.18
[40] lavaan_0.6-2          ellipsis_0.3.0        manipulate_1.0.1
[43] pkgconfig_2.0.1       nnet_7.3-12           deldir_0.1-15
[46] tidyselect_1.0.0      rlang_0.4.5           reshape2_1.4.3
[49] TeachingDemos_2.10    munsell_0.5.0         tools_3.5.1
[52] cli_1.0.0             generics_0.0.2        broom_0.5.0
[55] stringr_1.3.1         ModelMetrics_1.2.2.2  caTools_1.17.1.1
[58] purrr_0.3.3           mitml_0.3-6           nlme_3.1-137
[61] lavaan.survey_1.1.3.1 compiler_3.5.1        curl_3.2
[64] tibble_3.0.0          pbivnorm_0.6.0        stringi_1.2.4
[67] TAM_2.12-18           uroot_2.0-9           nloptr_1.0.4
[70] vegan_2.5-2           permute_0.9-4         urca_1.3-0
[73] vctrs_0.2.4           pillar_1.4.3          LearnBayes_2.15.1
[76] lifecycle_0.2.0       data.table_1.11.4     bitops_1.0-6
[79] R6_2.2.2              KernSmooth_2.23-15    codetools_0.2-15
[82] boot_1.3-20           gtools_3.8.1          assertthat_0.2.0
[85] withr_2.1.2           fracdiff_1.4-2        mnormt_1.5-5
[88] Deriv_3.8.5           mgcv_1.8-24           expm_0.999-2
[91] quadprog_1.5-5        grid_3.5.1            rpart_4.1-13
[94] tidyr_1.0.2           coda_0.19-1           class_7.3-14
[97] minqa_1.2.4           lubridate_1.7.4
>
