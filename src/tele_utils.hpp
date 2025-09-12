#ifndef TELE_UTILS_HPP
#define TELE_UTILS_HPP

#include <algorithm>
#include <cmath>
#include <vector>

#include <Rcpp.h>
using namespace Rcpp;

/*
 * Compute (exp(x * y) - 1) / x, preserving numerical accuracy for x close to
 * zero. Moreover, return y when x = 0, consistent with the asymptotic.
 */
inline double exp_ratio(double x, double y)
{
  return x != 0.0 ? std::expm1(x * y) / x : y;
}

/*
 * Let x be a vector containing the turning points of a telegraph process X(t)
 * from time 0 up to max(x). They need not be sorted, and 0 must not be
 * included.
 *
 * (i) If antithetic == false, compute
 *
 *        S0 * max(A - K, 0) ,
 *
 * where A is
 *
 *        exp(a * X(t) + b)
 *
 * integrated with respect to t from 0 to max(x).
 *
 * (ii) If antithetic == true, compute
 *
 *        0.5 * S0 * max(A - K, 0) + 0.5 * S0 * max(Aa - K, 0) ,
 *
 * where A and Aa are, respectively,
 *
 *        exp(a * X(t) + b)   and   exp(-a * X(t) + b)
 *
 * integrated with respect to t from 0 to max(x).
 */
template<bool antithetic>
double
  tele_payoff(const NumericVector& x, double a, double b, double S0, double K)
{
  int n = x.size();
  std::vector<double> xc(n);

  std::copy(x.begin(), x.end(), xc.begin());
  std::sort(xc.begin(), xc.end());

  double S = xc[0], A = exp_ratio(a + b, xc[0]), par = -1.0;

  if constexpr (antithetic) {
    double Aa = exp_ratio(-a + b, xc[0]);

    for (int i = 1; i < n; ++i) {
      A += exp_ratio(par * a + b, xc[i] - xc[i - 1])
           * std::exp(a * S + b * xc[i - 1]);
      Aa += exp_ratio(par * -a + b, xc[i] - xc[i - 1])
            * std::exp(-a * S + b * xc[i - 1]);
      S += par * (xc[i] - xc[i - 1]);
      par *= -1.0;
    }

    return 0.5 * std::max(S0 * A - K, 0.0) + 0.5 * std::max(S0 * Aa - K, 0.0);
  } else {
    for (int i = 1; i < n; ++i) {
      A += exp_ratio(par * a + b, xc[i] - xc[i - 1])
           * std::exp(a * S + b * xc[i - 1]);
      S += par * (xc[i] - xc[i - 1]);
      par *= -1.0;
    }

    return std::max(S0 * A - K, 0.0);
  }
}

#endif // TELE_UTILS_HPP