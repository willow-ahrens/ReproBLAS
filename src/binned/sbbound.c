#include <math.h>

#include <binned.h>

#include "../common/common.h"

/**
 * @brief Get binned single precision summation error bound
 *
 * This is a bound on the absolute error of a summation using binned types
 *
 * @param fold the fold of the binned types
 * @param N the number of single precision floating point summands
 * @param X the summand of maximum absolute value
 * @param S the value of the sum computed using binned types
 * @return error bound
 *
 * @author Willow Ahrens
 * @date   31 Jul 2015
 *
 * @author Willow Ahrens
 * @author Hong Diep Nguyen
 * @date   21 May 2015
 */
float binned_sbbound(const int fold, const int N, const float X, const float S) {
  return (float)(MAX(fabs((double)X), ldexp(0.5, FLT_MIN_EXP - 1)) * ldexp(0.5, (1 - fold) * SBWIDTH + 1) * N + ((7.0 * FLT_EPSILON) / (1.0 - 6.0 * sqrt((double)FLT_EPSILON) - 7.0 * FLT_EPSILON)) * fabs((double)S));
}
