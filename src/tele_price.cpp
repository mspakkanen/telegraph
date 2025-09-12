// [[Rcpp::plugins(cpp17)]]

#include <cmath>

#include <Rcpp.h>
#include "tele_utils.hpp"
using namespace Rcpp;

/*
 * Compute a Monte Carlo estimate of
 *
 *      E[max(A - K, 0)] ,
 *
 * where A is the integral with respect to t of
 *
 *      S0 * exp(sigma * X(t) - 0.5 * sigma^2)
 *
 * over [0, 1], where X(t) is a telegraph process with underlying Poisson
 * intensity lambda. Use N mutually independent paths.
 *
 * Reduce bias by:
 *
 * (i) multiplying X(t) with an independent Rademacher random variable (-1 with
 * prob 0.5 and 1 with prob 0.5) if antithetic == false,
 *
 * (ii) using antithetic pairs of paths if antithetic == true.
 */

// [[Rcpp::export]]
double tele_price(
  double S0, double K, double sigma, double lambda, int N, bool antithetic)
{
  IntegerVector x(N);
  x = rpois(N, lambda);

  NumericVector S(N), U(max(x) + 1);

  for (int i = 0; i < N; ++i) {
    if (x[i] > 0) {
      U[Range(0, x[i] - 1)] = runif(x[i], 0, 1);
    }
    U[x[i]] = 1;
    if (antithetic) {
      S[i] = tele_payoff<true>(U[Range(0, x[i])], sigma * std::sqrt(lambda),
                               -0.5 * sigma * sigma, S0, K);
    } else {
      S[i] = tele_payoff<false>(U[Range(0, x[i])],
                                (2 * rbinom(1, 1, 0.5)[0] - 1) * sigma
                                  * std::sqrt(lambda),
                                -0.5 * sigma * sigma, S0, K);
    }
  }

  return mean(S);
}