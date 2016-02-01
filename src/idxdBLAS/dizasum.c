#include "idxdBLAS.h"

/**
 * @brief Add to indexed double precision Y the absolute sum of complex double precision vector X
 *
 * Add to Y the indexed sum of magnitudes of elements of X.
 *
 * @param fold the fold of the indexed types
 * @param N vector length
 * @param X complex double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y indexed scalar Y
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void idxdBLAS_dizasum(const int fold, const int N, const void *X, const int incX, double_indexed *Y){
  idxdBLAS_dmzasum(fold, N, X, incX, Y, 1, Y + fold, 1);
}
