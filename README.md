# Bayesian Analysis of Penguin Bill Length

This repository contains the code and Stan models used for the Bayesian analysis of penguin bill length data using the Palmer Penguins dataset.

## Project Overview

The project compares six Bayesian regression models for predicting penguin bill length:

### Normal Models

1. Pooled Normal Model
2. Separate Normal Model
3. Hierarchical Normal Model

### Student-t Models

4. Pooled Student-t Model
5. Separate Student-t Model
6. Hierarchical Student-t Model

The models were implemented in Stan and fitted in R using the `rstan` package.

## Repository Structure

```text id="sk2v9m"
.
├── analysis.R
├── stan/
│   ├── normal_pooled.stan
│   ├── normal_separate.stan
│   ├── normal_hierarchical.stan
│   ├── student_t_pooled.stan
│   ├── student_t_separate.stan
│   └── student_t_hierarchical.stan
└── README.md
```

## Data

The analysis uses the Palmer Penguins dataset.
The response variable is:

* `bill_length_mm`

Predictors include:

* Standardized bill depth
* Species
* Sex

## Software and Packages

The project was developed in R.

Main packages used:

```r id="bz39ka"
library(rstan)
library(loo)
library(bayesplot)
library(posterior)
library(ggplot2)
library(dplyr)
```

## Model Description

### Pooled Models

Assume all penguin species share the same regression parameters.

### Separate Models

Estimate independent intercepts for each species.

### Hierarchical Models

Assigning hierarchical priors to species-level parameters.

### Student-t Models

Replace the normal likelihood with a Student-t likelihood to provide robustness against outliers.

## Bayesian Framework

Weakly informative priors were used for the regression parameters and variance terms.

Model comparison was conducted using:

* Posterior predictive checks
* Trace plots
* Convergence diagnostics
* Leave-One-Out Cross Validation (LOO-CV)

## Running the Analysis

1. Open the R script.
2. Make sure all required packages are installed.
3. Ensure the Stan files are located in the `stan/` folder.
4. Run the analysis script.

## Author

Sima Albalkhi
MSc Statistics and Data Science
Johannes Kepler University Linz

## Course Information

Bayesian Statistics
Summer Semester 2026
