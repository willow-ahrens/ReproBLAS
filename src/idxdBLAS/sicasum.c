/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "idxdBLAS.h"

/**
 * @brief Compute indexed single precision absolute sum Z of complex single precision vector X
 *
 * Set Z to the indexed sum of magnitudes of elements of X.
 *
 * @param N vector length
 * @param X complex single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Z indexed scalar Z
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void idxdBLAS_sicasum(const int fold, const int N, const void *X, const int incX, float_indexed *Y){
  idxdBLAS_smcasum(fold, N, X, incX, Y, 1, Y + fold, 1);
}
