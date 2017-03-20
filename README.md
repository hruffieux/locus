# locus – large-scale variational inference for combined predictor and outcome selection

[![Travis-CI Build Status](https://travis-ci.org/hruffieux/locus.svg?branch=master)](https://travis-ci.org/hruffieux/locus)
 
## Overview

**locus** is an R package providing efficient variational algorithms for
simultaneous variable selection of predictors and associated outcomes based
on multivariate regression models. Dependence across outcomes linked to the 
same predictors is captured through the model hierarchical structure 
(H. Ruffieux, A. C. Davison, J. Hager, I. Irincheeva, Efficient inference 
for genetic association studies with multiple outcomes, *Biostatistics*, 2017). 

## Installation

To install, run the following commands in R:

``` r
install.packages("devtools")
devtools::install_github("hruffieux/locus")
```

## Algorithms

The algorithms for joint predictor and response selection provided in **locus**
implement inference for regression models with 

* identity link;
* logistic link;
* probit link; 
* identity-probit link.

Moreover, predictor-level external knowledge can be incorportated to inform the
selection.

## License and authors

This software uses the GPL v2 license, see [LICENSE](LICENSE).
Authors and copyright are provided in [DESCRIPTION](DESCRIPTION). Loris Michel
has also contributed to the development of this project.

## Issues

To report an issue, please use the [locus issue tracker](https://github.com/hruffieux/locus/issues) at github.com.
