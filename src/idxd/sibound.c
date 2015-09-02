#include <math.h>

#include <idxd.h>

#include "../common/common.h"

/**
 * @brief Get indexed single precision summation error bound
 *
 * This is a bound on the absolute error of a summation using indexed types
 *
 * @param fold the fold of the indexed types
 * @param N the number of single precision floating point summands
 * @param X the summand of maximum absolute value
 * @param S the value of the sum computed using indexed types
 * @return error bound
 *
 * @author Peter Ahrens
 * @date   31 Jul 2015
 *
 * @author Peter Ahrens
 * @author Hong Diep Nguyen
 * @date   21 May 2015
 */
float idxd_sibound(const int fold, const int N, const float X, const float S) {
  return (float)(MAX(fabs((double)X), ldexp(0.5, FLT_MIN_EXP - 1)) * ldexp(0.5, (1 - fold) * SIWIDTH + 1) * N + ((7.0 * FLT_EPSILON) / (1.0 - 6.0 * sqrt((double)FLT_EPSILON) - 7.0 * FLT_EPSILON)) * fabs((double)S));
}
