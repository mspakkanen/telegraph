# Simulation methods for telegraph processes

This repo contains C++/Rcpp simulation code from the paper:

* G. Barrera, J. Lukkarinen and M. S. Pakkanen (2025): Wasserstein error estimates between telegraph processes and Brownian motion, [arXiv:2509.11871](https://arxiv.org/abs/2509.11871). 

## Contents and getting started

```
.
├── src
│   ├── bm_price.cpp   # Arithmetic Asian option Monte Carlo pricer under Brownian motion
│   ├── tele_price.cpp # Arithmetic Asian option Monte Carlo pricer under a telegraph process
│   └── tele_utils.hpp # Utility functions to compute the option payoff for a telegraph process
├── LICENSE            # MIT License
└── README.md          # This file
```

Copy the contents of `src/` to your working directory and `Rcpp::sourceCpp()` the files `bm_price.cpp` and `tele_price.cpp`.

## Requirements

* [R](https://www.r-project.org/)
* [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) and the requisite C++ compiler tool chain
