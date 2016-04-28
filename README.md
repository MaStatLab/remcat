# remcat - REcursive Mixtures for Conditional Association Testing

This package carries out an adaptive test for conditional association in contigency tables using recursive mixture models.

### Install
The package can be installed on Linux and Mac using `devtools`:

```S
install.packages('devtools')
library('devtools')
devtools::install_github('remcat', 'MaStatLab')
```

### Use

```S
remcat = function( xg.mat, xe.mat, y.vec, rho0)
```

### Reference

Ma L. (2013) Adaptive Testing of Conditional Association Through Recursive Mixture Modeling. Journal of the American Statistical Association. Volume 108, Issue 504, pages 1493-1505.
