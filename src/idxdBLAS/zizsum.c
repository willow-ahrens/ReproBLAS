#include "idxdBLAS.h"

/**
 * @brief Add to indexed complex double precision Y the sum of complex double precision vector X
 *
 * Add to Y the indexed sum of X.
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
void idxdBLAS_zizsum(const int fold, const int N, const void *X, const int incX, double_complex_indexed *Y){
  idxdBLAS_zmzsum(fold, N, X, incX, Y, 1, Y + 2 * fold, 1);
}
