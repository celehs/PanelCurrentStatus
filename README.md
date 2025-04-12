# PanelCurrentStatus: Risk Prediction Models with Panel Current Status Data

## Overview

The `PanelCurrentStatus` package implements methods for developing and evaluating risk prediction models using panel current status data. In panel current status data, subjects are examined for the occurrence of an event at several pre-scheduled visit times rather than being continuously monitored. Such data arise frequently in biomedical studies where continuous monitoring is impractical or costly.

This package provides tools to:
- Compute the conditional censoring logistic (CCL) estimator, which transforms panel current status data into a binary outcome analysis
- Evaluate prediction performance of estimated risk models with panel current status data
- Calculate model metrics from ROC curves via kernel smoothing

The CCL estimator offers advantages over existing methods by:
- Building on logistic regression estimators while incorporating monitoring time information
- Providing better performance in relatively small sample sizes
- Remaining robust under various model specifications

## Installation

Install development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("celehs/PanelCurrentStatus")
```

## Key Functions

- `ccl.fit()`: Computes the scaled coefficients from the conditional censoring logistic estimator
- `ccl.roc()`: Evaluates model metrics from ROC curve via kernel smoothing

## References

Chan S, Wang X, Jazić I, Peskoe S, Zheng Y, Cai T. Developing and evaluating risk prediction models with panel current status data. Biometrics. 2021 Jun;77(2):599-609. doi: 10.1111/biom.13317. Epub 2020 Jul 8. PMID: 32562264; PMCID: PMC8168594.
