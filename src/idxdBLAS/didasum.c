/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "idxdBLAS.h"

/**
 * @brief Add to indexed double precision Y the absolute sum of double precision vector X
 *
 * Add to Y the indexed sum of absolute values of elements in X.
 *
 * @param N vector length
 * @param X double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y indexed scalar Y
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void idxdBLAS_didasum(const int fold, const int N, const double *X, const int incX, double_indexed *Y){
  idxdBLAS_dmdasum(fold, N, X, incX, Y, 1, Y + fold, 1);
}
