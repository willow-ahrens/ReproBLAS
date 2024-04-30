#include <binned.h>

/**
 * @brief  Add binned single precision vectors (Y += X)
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
void binned_sbsbaddv(const int fold, const int N, const float_binned *X, const int incX, float_binned *Y, const int incY){
  int i;
  for(i = 0; i < N; i++, X += incX * binned_sbnum(fold), Y += incY * binned_sbnum(fold)){
    binned_smsmadd(fold, X, 1, X + fold, 1, Y, 1, Y + fold, 1);
  }
}
