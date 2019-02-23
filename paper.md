---
title: 'ChiRP: Chinese Restaurant Process Mixtures for Regression and Clustering'
tags:
  - R
  - Bayesian
  - Nonparametric
  - Clustering
  - Dirichlet Process
  - Chinese Restaurant
authors:
  - name: Arman Oganisian
    orcid: 0000-0002-0437-4611
    affiliation: 1
affiliations:
 - name: Department of Biostatistics, Epidemiology, and Informatics, University of Pennsylvania
   index: 1
date: 15 February 2019
bibliography: paper.bib
nocite: | 
  @hannah2011, @muller2015, @Neal2000, @gershman2012, @stephens2000, @Rod2014, @Hastings1970
---

# Summary and Innovation
Regression modeling and clustering are common statistical tasks in biomedical research. However, regression problems often involve parametric assumptions (e.g. Normality). Clustering often involves pre-specifying the number of clusters - a quantity typically unknown to the researcher. Flexible machine learning methods do exist for such problems, but they tend to focus on point estimation rather than interval estimation. The latter is often just as important as the former in biomedical settings. There is great need for regression methods that do not make rigid parametric assumptions and clustering methods that do not require pre-specifying the number of clusters. It would also be ideal for these methods to produce uncertainty estimates that allow for valid statistical inference in addition to predictions.

`ChiRP` package implements **Chi**nese **R**estaurant **P**rocess (CRP) mixtures in `R`. CRP mixtures are a flexible class of Bayesian nonparametric models that can be used for general regression modeling and clustering problems that meets the needs mentioned above. 

At a high level, the model works by partitioning a complex, heterogenous dataset into more homogenous clusters - each associated with a locally parametric regression model. Unlike traditional clustering procedures, CRP mixtures allow for infinitely many clusters - thus removing the need to specify the number of clusters in advance. The CRP clustering also allows us to train nonparametric regressions that generate testing and training predictions by ensembling over the local cluster-specific regression models. Probabilistically valid interval estimates are also given for each prediction. We note that while the local cluster-specific models are parametric, allowing for infinitely many clusters yields an infinite-dimensional parameter space.

# Outcome Types and Model Output

We may be given a training dataset with $n$ subjects $D_T=(Y_i, X_i)_{i=1:n}$. Here, $Y_i$ is the scalar outcome/label and $X_i$ is a $p\times1$ feature vector.`ChiRP` trains a fully Bayesian CRP model with Monte Carlo Markov Chain (MCMC) and yields the following output:

1. In-sample predictions $( \hat{Y}_{i} )_{i=1:n}$ from a nonparametric CPR regression of $Y$ on $X$.
2. Out-of-sample predictions on an unlabeled test set $(\tilde{X}_i)_{1:m}$, $(\hat{\tilde{Y}}_i)_{i=1:m}$.
3. Latel cluster membership for both training and testing subjects.

`ChiRP`'s implementation is fully Bayesian - yielding not only point predictions, but interval estimtes that allow us to measure uncertainty around these predictions. Uncertainty estimation is often just as important, if not more important, than point estimation in biomedical sciences. 

`ChiRP` implements three different local cluster-specific regressions depending on the nature of $Y_i$: 

1. Continuous outcomes using cluster-specific linear regressions.
2. Binary outcomes using cluster-specific logistic regressions.
3. Zero-inflated, semi-continuous outcomes using cluster-specific zero-inflated regressions. [See @Oganisian2018]



# References
