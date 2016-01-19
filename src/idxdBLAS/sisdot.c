#include "idxdBLAS.h"

/**
 * @brief Add to indexed single precision Z the dot product of single precision vectors X and Y
 *
 * Add to Z the indexed sum of the pairwise products of X and Y.
 *
 * @param fold the fold of the indexed types
 * @param N vector length
 * @param X single precision vector
 * @param incX X vector stride (use every incX'th element)
 * @param Y single precision vector
 * @param incY Y vector stride (use every incY'th element)
 * @param Z indexed scalar Z
 *
 * @author Peter Ahrens
 * @date   15 Jan 2016
 */
void idxdBLAS_sisdot(const int fold, const int N, const float *X, const int incX, const float *Y, const int incY, float_indexed *Z){
  idxdBLAS_smsdot(fold, N, X, incX, Y, incY, Z, 1, Z + fold, 1);
}
