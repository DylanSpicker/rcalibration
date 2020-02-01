# rcalibration
Implements Generalized Regression Calibration (and Related) Estimators for Measurement Error Correction.

This package implements a generalized regression calibration estimator, for the correction of measurement error in standard modelling contexts. This generalized estimator contains the more common replicate-based estimators (Carroll, et al., 2006) as a special case and incorporates a broad correction that establishes consistency under a broad class of error models. The package further provides utilities related to these corrections, notably the capacity to compute optimal weights when taking a convex combination of error-prone proxies and utilities for estimating the structure of errors observed in data. 

## Installation
Install the latest version from github. Note, this requires `devtools`. 
```{r}
install.packages("devtools")
devtools::install_github("dylanspicker/rcalibration")
```

## Usage
TBD.