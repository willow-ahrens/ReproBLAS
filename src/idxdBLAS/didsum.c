/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "idxdBLAS.h"

/**
 * @brief Compute indexed double precision sum Z of double precision vector X
 *
 * Set Z to the indexed sum of X.
 *
 * @param N vector length
 * @param X double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Z indexed scalar Z
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void idxdBLAS_didsum(const int fold, const int N, const double *X, const int incX, double_indexed *Y){
  idxdBLAS_dmdsum(fold, N, X, incX, Y, 1, Y + fold, 1);
}
