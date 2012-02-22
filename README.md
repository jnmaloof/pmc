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





