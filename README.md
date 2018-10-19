## LOCUS: Large-scale variational inference for Bayesian variable selection in multiple-response regression

[![Travis-CI Build Status](https://travis-ci.org/hruffieux/locus.svg?branch=master)](https://travis-ci.org/hruffieux/locus)

## Overview

**LOCUS** is an R package providing efficient variational algorithms for
simultaneous selection of responses and associated predictors based
on hierarchical sparse regression models. The method can for instance be 
used for genome-wide association analyses with multiple clinical endpoints
or for molecular quantitative trait locus (QTL) analyses which involve 
hundreds of thousands of genetic variants as candidate predictors and
thousands of molecular levels as responses, for thousands of samples.  

**LOCUS** is described and extensively assessed in H. Ruffieux, 
A. C. Davison, J. Hager and I. Irincheeva, Efficient inference
for genetic association studies with multiple outcomes, *Biostatistics*, 
18:618:636, 2017, doi: 10.1093/biostatistics/kxx007.

## Installation

To install, run the following commands in R:

``` r
require(devtools) # after having installed devtools (install.packages("devtools"))
devtools::install_github("hruffieux/locus")
```
## License and authors

This software uses the GPL v2 license, see [LICENSE](LICENSE).
Authors and copyright are provided in [DESCRIPTION](DESCRIPTION). Loris Michel
has also contributed to the development of this project.

## Issues

To report an issue, please use the [locus issue tracker](https://github.com/hruffieux/locus/issues) at github.com.