#include <math.h>

#include <indexed.h>

/**
 * @brief Get indexed double precision summation error bound
 *
 * This is a bound on the absolute error of a summation using indexed types
 *
 * @param fold the fold of the indexed types
 * @param N the number of double precision floating point summands
 * @param X the maximum absolute value of the summands
 * @return error bound
 *
 * @author Peter Ahrens
 * @author Hong Diep Nguyen
 * @date   21 May 2015
 */
double dibound(const int fold, const int N, const double X) {
  return X * ldexp(0.5, (1 - fold)*(DIWIDTH - 1) + 1) * N;
}
