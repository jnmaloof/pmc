## Phylogenetic Monte Carlo

* Author: Carl Boettiger
* License: [CC0](http://creativecommons.org/publicdomain/zero/1.0/)

This package accompanies my publication, with Graham Coop & Peter Ralph (2012) **Is your phylogeny 
Informative? Measuring the power of comparative methods**. In *Evolution*.
([doi](http://dx.doi.org/10.1111/j.1558-5646.2012.01574.x)) 
([pdf](http://www.mendeley.com/download/public/98752/4485545653/9a209c7dd29980fd2e47c06eb8b2d1d7dd6f70d4/dl.pdf))
([arXiv](http://arxiv.org/abs/1110.4944))
([code](https://github.com/cboettig/pmc))
([data](http://datadryad.org/handle/10255/dryad.37645))



## Install

The easiest way to install the current release is through the CRAN repository, 

```r
install.packages("pmc")
```

To install the current development version, there are a few options.  One is to install `devtools` package, which will let you install directly from github:

```r
install.packages("devtools")
library(devtools)
install_github("pmc", "cboettig")
```

The other option is to download the current tarball, put it in the working directory, and then install from R using
```r
install.packages("pmc*.tar.gz", ,type="source", repos=NULL)
```

## Getting Started
The best place to start learning the package commands is to read the [package vignette](https://github.com/cboettig/pmc/blob/master/vignettes/pmc_tutorial.pdf), copy-pasting the examples in as you go.  For a faster start, just look at the examples in the function documentation, i.e. `?pmc`.  

## Bug Reports
If you find the package is giving unexpected results, is missing a feature you'd like, or have any other questions about the package, don't hesitate to email me. 

You can also [file an issue](https://github.com/cboettig/pmc/issues) which will provide a record of the bug/request.  Posting an issue also automatically sends me an email.  





