
# PanelCurrentStatus: Risk Prediction Models with Panel Current Status Data

[![CRAN](https://www.r-pkg.org/badges/version/PanelCurrentStatus)](https://CRAN.R-project.org/package=PanelCurrentStatus)

## Overview

This package contains R functions to compute the conditional censoring
logistic (CCL) estimator and model metrics to evaluate risk predictions
using panel current status data. The CCL estimator takes advantage of
the ability to transform panel current status data into a binary outcome
analysis, building on existing logistic regression estimators by
incorporating monitoring time information into the working model.

## Installation

Install development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("celehs/PanelCurrentStatus")
```

## Getting Started

Please click
[HERE](https://github.com/celehs/PanelCurrentStatus/blob/master/demo.pdf)
to view a demo. The data example in the demo can be downloaded
[HERE](https://github.com/celehs/PanelCurrentStatus/blob/master/data.rds).

## References

Chan S, Wang X, Jazić I, Peskoe S, Zheng Y, Cai T. Developing and
evaluating risk prediction models with panel current status data.
Biometrics. 2020 Jun 19. doi: 10.1111/biom.13317. Epub ahead of print.
PMID: 32562264.
