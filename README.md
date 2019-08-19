# Risk Prediction Models with Panel Current Status Data

## Overview

This package contains R functions to compute the conditional censoring logistic (CCL) estimator and model metrics to evaluate risk predictions using panel current status data. The CCL estimator takes advantage of the ability to transform panel current status data into a binary outcome analysis, building on existing logistic regression estimators by incorporating monitoring time information into the working model. 

## Installation

Install package from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("celehs/PanelCurrentStatus")
```

## Getting Started

Please click [HERE](https://github.com/celehs/PanelCurrentStatus/blob/master/demo.pdf) to view a demo.  The data example in the demo can be downloaded [HERE](https://github.com/celehs/PanelCurrentStatus/blob/master/data.rds).

## References

S. Chan, X. Wang, I. Jazic, S. Peskoe, Y. Zheng, T. Cai. Developing and Evaluating Risk Prediction Models with Panel Current Status Data. Submitted to _Biometrics_.
