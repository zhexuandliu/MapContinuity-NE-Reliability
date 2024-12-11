
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Modified Rtsne() function that returns the similarity matrix in input space

## Installation

This package is modified from R package `Rtsne` (version 0.17) \[1\]. We
added an extra output to have access to the similarity matrix in the
input space, i.e. $$p_{ij}$$’s \[2\].

To install the package from the github repository, use:

``` r
if(!require(devtools)) install.packages("devtools") # If not already installed
devtools::install_github("zhexuandliu/RtsneWithP")
```

For Mac users, it is possible to encounter error when installing the package `RtsneWithP` if `gfortran` is not already installed. Please see https://github.com/zhexuandliu/MapContinuity-NE-Reliability/blob/main/RtsneWithP/gfortran-install.md for detailed instructions on the installation of `gfortran`. Or install the package from binaries, which can be downloaded from https://github.com/zhexuandliu/MapContinuity-NE-Reliability/tree/main/Code%20and%20data/PackageBinaries.

## Usage

``` r
library(RtsneWithP) # Load package
iris_unique <- unique(iris) # Remove duplicates
set.seed(42) # Sets seed for reproducibility
tsne_out <- RtsneWithP::Rtsne(as.matrix(iris_unique[,1:4]), theta = 0) # set theta=0 to run exact tSNE
print(tsne_out$P[1:5,1:5])
#>               [,1]          [,2]          [,3]          [,4]          [,5]
#> [1,] 5.173550e-312  9.085358e-05  9.642498e-05  5.262639e-05  4.233656e-04
#> [2,]  9.085358e-05 5.259419e-312  3.056671e-04  2.825609e-04  6.190919e-05
#> [3,]  9.642498e-05  3.056671e-04 5.537960e-312  3.520082e-04  1.020901e-04
#> [4,]  5.262639e-05  2.825609e-04  3.520082e-04 4.928655e-312  5.495371e-05
#> [5,]  4.233656e-04  6.190919e-05  1.020901e-04  5.495371e-05 5.638258e-312
```

## Details

This R package offers an extra return value based on the Rtsne()
function of the original package \[1\]. Note that the returned
similarity matrix is accurate only when the parameter `theta` is set to
be $0$ (exact tSNE).

## References

\[1\] <https://github.com/jkrijthe/Rtsne>

\[2\] L.J.P. van der Maaten and G.E. Hinton. “Visualizing
High-Dimensional Data Using t-SNE.” Journal of Machine Learning Research
9(Nov):2579-2605, 2008.
