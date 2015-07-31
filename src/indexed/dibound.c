#include <math.h>

#include <indexed.h>

#include "../common/common.h"

/**
 * @brief Get indexed double precision summation error bound
 *
 * This is a bound on the absolute error of a summation using indexed types
 *
 * @param fold the fold of the indexed types
 * @param N the number of double precision floating point summands
 * @param X the summand of maximum absolute value
 * @param S the value of the sum computed using indexed types
 * @return error bound
 *
 * @author Peter Ahrens
 * @date   31 Jul 2015
 */
double dibound(const int fold, const int N, const double X, const double S) {
  double g = ((7.0 - 2.0 * sqrt(DBL_EPSILON)) / ((1.0 - 2.0 * sqrt(DBL_EPSILON)) * (1.0 - 4.0 * sqrt(DBL_EPSILON)))) * DBL_EPSILON;
  return MAX(fabs(X), ldexp(0.5, DBL_MIN_EXP - 1)) * ldexp(0.5, (1 - fold)*(DIWIDTH - 1) + 1) * N + (g / (1.0 - g)) * fabs(S);
}
