# Repository for `ChiRP` Package <a href="url"><img src="logo.png" align="right" height="200" width="200" ></a>


## About
The R package `ChiRP` is an MCMC-based implementation of **Chi**nese **R**estaurant **P**rocess (CRP) Mixtures of various regression and clustering analysis. 

CRP models (aka Dirichlet Process models) are a class of Bayesian nonparametric models.

Install using `devtools` package
```
## install.packages('devtools' ) ## make sure to have devtools intalled 
devtools::install_github('stablemarkets/ChiRP')
library(ChiRP)
``` 

Note this package is still in beta. You can follow updates about the package on [twitter](https://twitter.com/StableMarkets) .

## Companion Web Site, Examples, and Documentation
The companion web site at (https://stablemarkets.github.io/ChiRPsite/index.html) contains the statistical details of the model as well as several examples. 

Help documentation in `R` is also available.  After installing the package and loading it with `library()`, use `?` to access help documentation for specific functions:
```
?ChiRP::NDPMix
?ChiRP::ZDPMix
?ChiRP::PDPMix
?ChiRP::cluster_assign_mode
``` 
The help file for each section contains an example that you can copy/paste and run directly in your `R` session.

## How to Cite this Package
Since this package was written in conjunction with [this paper](https://arxiv.org/abs/1810.09494), please cite it when using this package. Thank you!

## Reporting Bugs and Issues
If you have any issues with the code, please open an issue on GitHub.

## Contact
The corresponding package author is Arman Oganisian (email: aoganisi@upenn.edu).