# Simulation methods for telegraph processes

This repo contains C++/Rcpp simulation code from the paper:

* G. Barrera, J. Lukkarinen and M. S. Pakkanen (2025): Wasserstein error estimates between telegraph processes and Brownian motion.

## Contents

```
.
├── src
|   ├── bm_price.cpp   # Arithmetic Asian option Monte Carlo pricer under Brownian motion
|   ├── tele_price.cpp # Arithmetic Asian option Monte Carlo pricer under a telegraph process
|   └── tele_utils.hpp # Utility functions to compute the option payoff for a telegraph process
├── LICENCE
└── README.md          # This file
```

## Requirements

* [R](https://www.r-project.org/)
* [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) and the requisite C++ compiler tool chain
