#include "idxdBLAS.h"

/**
 * @brief Add to indexed complex single precision Y the sum of complex single precision vector X
 *
 * Add to Y the indexed sum of X.
 *
 * @param N vector length
 * @param X complex single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param indexed scalar Y
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void idxdBLAS_cicsum(const int fold, const int N, const void *X, const int incX, float_complex_indexed *Y){
  idxdBLAS_cmcsum(fold, N, X, incX, Y, 1, Y + 2 * fold, 1);
}
