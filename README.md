# Predicting eczema severity using serum biomarkers

This repository contains the code for the article by [**Hurault et al. (preprint), "Can serum biomarkers predict the outcome of systemic therapy for atopic dermatitis?"**]().

The code is written in the R language for statistical computing and the models using the probabilistic programming language [Stan](https://mc-stan.org/).

## File structure

The dataset used in this study is not available according to our data sharing agreement.
During the analysis, the dataset is loaded from a proprietary package `TanakaData` which includes the raw files as well as data processing functions.

Utility functions used within the scripts are available in [`functions.R`](functions.R).
In addition, we used functions from Guillem Hurault's personal package, [HuraultMisc](https://github.com/ghurault/HuraultMisc).

The `Models` folder contains the different Stan models developed in this project:

- [`RW.stan`](Models/RW.stan): the random walk model, one of the reference model.
- [`AR.stan`](Models/AR.stan): the autoregressive model, one of the reference model.
- [`MixedAR.stan`](Models/MixedAR.stan): the mixed effect autoregressive model, one of the reference model.
- [`SSM.stan`](Models/SSM.stan): the Bayesian state space model without covariates.
- [`SSMX.stan`](Models/SSMX.stan): the Bayesian state space model with covariates (following a horseshoe prior).

The modelling workflow in separated into different scripts:

- [`check_models.R`](check_models.R): Conduct prior predictive checks and fake data check of the different models.
This script is notably useful to simulate data that resembles the one we used.
- [`fit_models.R`](fit_models.R): Fit the different models to real data, perform diagnostics and posterior predictive checks.
- [`run_validation.R`](run_validation.R): Run the validation process (K-fold cross-validation and forward chaining).
- [`check_performance.R`](check_performance.R): Analyse validation results.

## License

This open source version of this project is licensed under the GPLv3 license, which can be seen in the [LICENSE](LICENSE) file.
