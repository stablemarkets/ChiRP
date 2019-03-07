## `ChiRP`: Chinese Restaurant Process Mixtures for Regression and Clustering <a href="url"><img src="logo.png" align="right" height="200" width="200" ></a>

Development Status:

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.org/stablemarkets/ChiRP.svg?branch=master)](https://travis-ci.org/stablemarkets/ChiRP)
[![Coveralls github](https://img.shields.io/coveralls/github/stablemarkets/ChiRP.svg?style=popout)](https://coveralls.io/github/stablemarkets/ChiRP)

[![status](http://joss.theoj.org/papers/3b83a0a3f1220f97657a1075b78e480a/status.svg)](http://joss.theoj.org/papers/3b83a0a3f1220f97657a1075b78e480a)


## About
The R package `ChiRP` is an MCMC-based implementation of **Chi**nese **R**estaurant **P**rocess (CRP) mixtures for regression and clustering. CRP models (aka Dirichlet Process models) are a class of Bayesian nonparametric models. We provide facilities for zero-inflated semi-continuous outcomes, continuous outcomes, and binary outcomes.

## Installation

Install using `devtools` package
```
## install.packages('devtools' ) ## make sure to have devtools installed 
devtools::install_github('stablemarkets/ChiRP')
library(ChiRP)
``` 

## Documentation and Examples
The [companion web site](https://stablemarkets.github.io/ChiRPsite/index.html) contains the [statistical details](https://stablemarkets.github.io/ChiRPsite/modeldesc.html) of the model as well as several [replicable examples](https://stablemarkets.github.io/ChiRPsite/examples.html). 

Help documentation in `R` is also available.  After installing the package and loading it with `library()`, use `?` to access help documentation for specific functions:
```
?ChiRP::NDPMix  # for continuous outcomes
?ChiRP::ZDPMix  # for zero-inflated, semi-continuous outcomes
?ChiRP::PDPMix  # for binary outcomes
?ChiRP::cluster_assign_mode # computes posterior mode cluster assignment
``` 
The help file for each function contains an example that you can run directly in your `R` session.

## Reporting Issues
`ChiRP` uses the `testthat` package for unit-testing and Travis CI for continuous integration. Coverage of unit test is tracked using Coveralls. 

If you encounter any bugs or have feature requests, please [open an issue](https://github.com/stablemarkets/ChiRP/issues) on GitHub.

## Contributing to `ChiRP`
You can contribute in two ways:

1. Contribute to base code: First, start an issue in this repository with the proposed modification. Fork this repository, make changes/enhancements, then submit a pull request. The issue will be closed once the pull request is merged.
2. Contribute an example: First, start an issue in the [companion site's repository](https://github.com/stablemarkets/ChiRPsite). Fork the repository and add a new example to `examples.Rmd`. Use `rmarkdown::render_site()` to build the site. Submit a pull request in that same repository. The issue will be closed once updates are merged.

## Contact
The corresponding package author is Arman Oganisian (email: aoganisi@upenn.edu). You can follow updates about the package on [twitter](https://twitter.com/StableMarkets).
