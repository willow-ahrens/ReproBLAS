/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "idxdBLAS.h"

/**
 * @brief Compute indexed double precision absolute sum Z of complex double precision vector X
 *
 * Set Z to the indexed sum of magnitudes of elements of X.
 *
 * @param N vector length
 * @param X complex double precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Z indexed scalar Z
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void idxdBLAS_dizasum(const int fold, const int N, const void *X, const int incX, double_indexed *Y){
  idxdBLAS_dmzasum(fold, N, X, incX, Y, 1, Y + fold, 1);
}
