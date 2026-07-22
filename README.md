# Spatio-Temporal Covariance Model Fitting

This repository contains R scripts for estimating the parameters of a **spatio-temporal covariance model** 
from an empirical covariogram, using constrained nonlinear optimization (`nloptr`).
The spatio-temporal covariance model is built through the difference of two correlation 
models, characterized by negative values.

The estimation follows a **three-step fitting strategy**: the spatial and temporal 
marginal behaviors of the covariance are fitted first, and their estimates are then 
used as starting values (and constraints) for the joint spatio-temporal model. 
This sequential approach improves the stability and interpretability of the final 
optimization, which is typically harder to fit directly due to the larger parameter 
space and the nonlinear inequality constraints linking spatial, temporal, and 
interaction parameters.

## Repository contents

| File | Description |
|---|---|
| `Step_opt_params_S_T_ST.R` | Main script. Loads the empirical spatio-temporal covariance (cova_st.RData), produces 2D/3D plots, sources the three fitting steps in sequence, computes the fitted spatio-temporal covariance and plots the corresponding 3D surface, then computes goodness of fit (MAE, RMSE) indexes. |
| `Step1_C_space.R` | **Step 1 — Spatial marginal fitting.** Fits the spatial marginal covariance model to estimate parameters `A1, B1, c1, c2, c11, c21`, subject to specific constraints. |
| `Step2_C_time.R` | **Step 2 — Temporal marginal fitting.** Fits the temporal marginal covariance model to estimate parameters `A1, B1, a1, a2`, using `A1`/`B1` inherited from Step 1 and subject to specific constraints. |
| `Step3_C_spacetime.R` | **Step 3 — Joint spatio-temporal fitting.** Fits the spatio-temporal covariance model to estimate all parameters (`A1, B1, c1, c2, c11, c21, a1, a2`) jointly, using the Step 1 and Step 2 estimates as starting values, subject to specific constraints.


## Requirements

The scripts require the following R packages:

```r
install.packages(c("dplyr", "lattice", "RColorBrewer", "readxl", "nloptr"))
```

## Usage

1. Load the empirical spatio-temporal covariance data in the file `cova_st.RData`.
2. Run `Step_opt_params_S_T_ST.R`. This is the only script that needs to be executed 
directly — it sources `Step1_C_space.R`, `Step2_C_time.R`, and `Step3_C_spacetime.R` 
in the correct order.
3. Inspect the resulting plots and the `MAE`/`RMSE` values to assess model fit.

## Authors

- Prof. Sandra De Iaco — University of Salento, Italy
- Prof. Donato Posa — University of Salento, Italy
- Prof. Claudia Cappello — University of Salento, Italy
