#include <math.h>

#include <binned.h>

#include "../common/common.h"

/**
 * @brief Get binned double precision summation error bound
 *
 * This is a bound on the absolute error of a summation using binned types
 *
 * @param fold the fold of the binned types
 * @param N the number of double precision floating point summands
 * @param X the summand of maximum absolute value
 * @param S the value of the sum computed using binned types
 * @return error bound
 *
 * @author Peter Ahrens
 * @date   31 Jul 2015
 */
double binned_dbbound(const int fold, const int N, const double X, const double S) {
  return MAX(fabs(X), ldexp(0.5, DBL_MIN_EXP - 1)) * ldexp(0.5, (1 - fold) * DBWIDTH + 1) * N + ((7.0 * DBL_EPSILON) / (1.0 - 6.0 * sqrt(DBL_EPSILON) - 7.0 * DBL_EPSILON)) * fabs(S);
}
