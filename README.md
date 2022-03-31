# cdcatR: Cognitive Diagnostic Computerized Adaptive Testing in R
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/cdcatR?color=brightgreen)](https://cran.r-project.org/package=cdcatR)
[![](https://cranlogs.r-pkg.org/badges/cdcatR?color=blue)](https://cran.r-project.org/package=cdcatR)
[![](http://cranlogs.r-pkg.org/badges/grand-total/cdcatR?color=blue)](https://cran.r-project.org/package=cdcatR)

## How to cite the package

Sorrel, M.A., Nájera, P., and Abad, F.J. (2021). cdcatR: Cognitive Diagnostic Computerized Adaptive Testing in R. R package version 1.0.3 Retrieved from https://cran.r-project.org/package=cdcatR.

## Features of the package

* Conducting CD-CAT applications in R
* Comparing different CD-CAT applications in terms of classification accuracy and CAT length

## Installation

A stable version of `cdcatR` is available at CRAN and can be installed using:

```r
install.packages("cdcatR")
```

To install this package from source:

1. Windows users may need to install the [Rtools](https://cran.r-project.org/bin/windows/Rtools/) and include the checkbox option of installing Rtools to their path for easier command line usage. Mac users will have to download the necessary tools from the [Xcode](https://apps.apple.com/ca/app/xcode/id497799835?mt=12) and its related command line tools (found within Xcode's Preference Pane under Downloads/Components); most Linux distributions should already have up to date compilers (or if not they can be updated easily).
2. Install the `remotes` package (if necessary), and install the package from the Github source code.

```
#install.packages("remotes")
remotes::install_github("miguel-sorrel/cdcatR")
```