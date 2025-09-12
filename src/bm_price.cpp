#include <algorithm>
#include <cmath>

#include <Rcpp.h>
using namespace Rcpp;

/*
 * Compute a Monte Carlo estimate of
 *
 *      E[max(A - K, 0)] ,
 *
 * where A is the integral with respect to t of
 *
 *      S0 * exp(sigma * W(t) - 0.5 * sigma^2)
 *
 * on [0, 1], where W(t) is a standard Brownian motion. Discretise the
 * integral along M time steps. Use N mutually independent paths.
 *
 * Use antithetic pairs of paths if antithetic = TRUE.
 */

// [[Rcpp::export]]
double
  bm_price(double S0, double K, double sigma, int M, int N, bool antithetic)
{
  NumericVector drift(M + 1), X(M + 1), S(N);

  drift = Range(0, M);
  drift = -0.5 * sigma * sigma * drift / static_cast<double>(M);

  for (int i = 0; i < N; ++i) {
    X[Range(1, M)]
      = rnorm(M, 0.0, sigma * std::sqrt(1 / static_cast<double>(M)));
    X[0] = 0.0;
    std::partial_sum(X.begin(), X.end(), X.begin());

    if (antithetic) {
      S[i] = 0.5 * std::max(mean(S0 * exp(drift + X)) - K, 0.0)
             + 0.5 * std::max(mean(S0 * exp(drift - X)) - K, 0.0);
    } else {
      S[i] = std::max(mean(S0 * exp(drift + X)) - K, 0.0);
    }
  }

  return mean(S);
}