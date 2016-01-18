#include "idxdBLAS.h"

/**
 * @brief Add to indexed single precision Y the absolute sum of complex single precision vector X
 *
 * Add to Y the indexed sum of magnitudes of elements of X.
 *
 * @param N vector length
 * @param X complex single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y indexed scalar Y
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void idxdBLAS_sicasum(const int fold, const int N, const void *X, const int incX, float_indexed *Y){
  idxdBLAS_smcasum(fold, N, X, incX, Y, 1, Y + fold, 1);
}
