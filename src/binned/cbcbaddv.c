#include <binned.h>

/**
 * @brief  Add binned complex single precision vectors (Y += X)
 *
 * Performs the operation Y += X
 *
 * @param fold the fold of the binned types
 * @param N vector length
 * @param X binned vector X
 * @param incX X vector stride (use every incX'th element)
 * @param Y binned vector Y
 * @param incY Y vector stride (use every incY'th element)
 *
 * @author Willow Ahrens
 * @date   25 Jun 2015
 */
void binned_cbcbaddv(const int fold, const int N, const float_complex_binned *X, const int incX, float_complex_binned *Y, const int incY){
  int i;
  for(i = 0; i < N; i++, X += incX * binned_cbnum(fold), Y += incY * binned_cbnum(fold)){
    binned_cmcmadd(fold, X, 1, X + 2 * fold, 1, Y, 1, Y + 2 * fold, 1);
  }
}
