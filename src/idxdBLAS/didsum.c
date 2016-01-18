#include "idxdBLAS.h"

/**
 * @brief Add to indexed double precision Y the sum of double precision vector X
 *
 * Add to Y the indexed sum of X.
 *
 * @param N vector length
 * @param X double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y indexed scalar Y
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void idxdBLAS_didsum(const int fold, const int N, const double *X, const int incX, double_indexed *Y){
  idxdBLAS_dmdsum(fold, N, X, incX, Y, 1, Y + fold, 1);
}
