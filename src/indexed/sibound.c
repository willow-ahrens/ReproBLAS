#include <math.h>

#include <indexed.h>

/**
 * @brief Get indexed single precision summation error bound
 *
 * This is a bound on the absolute error of a summation using indexed types
 *
 * @param fold the fold of the indexed types
 * @param N the number of single precision floating point summands
 * @param X the maximum absolute value of the summands
 * @return error bound
 *
 * @author Peter Ahrens
 * @author Hong Diep Nguyen
 * @date   21 May 2015
 */
float sibound(const int fold, const int N, const float X) {
  return X * ldexpf(0.5, (1 - fold)*(SIWIDTH - 1) + 1) * N;
}
