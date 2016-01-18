#include "idxdBLAS.h"

/**
 * @brief Add to indexed complex single precision Z the conjugated dot product of complex single precision vectors X and Y
 *
 * Add to Z the indexed sum of the pairwise products of X and conjugated Y.
 *
 * @param N vector length
 * @param X complex single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y complex single precision vector
 * @param incY Y vector stride (use every incX'th element)
 * @param indexed scalar Z
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void idxdBLAS_cicdotc(const int fold, const int N, const void *X, const int incX, const void *Y, const int incY, float_complex_indexed *Z){
  idxdBLAS_cmcdotc(fold, N, X, incX, Y, incY, Z, 1, Z + 2 * fold, 1);
}
