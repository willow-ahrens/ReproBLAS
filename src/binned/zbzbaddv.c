#include <binned.h>

/**
 * @brief  Add binned complex double precision vectors (Y += X)
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
 * @author Peter Ahrens
 * @date   25 Jun 2015
 */
void binned_zbzbaddv(const int fold, const int N, const double_complex_binned *X, const int incX, double_complex_binned *Y, const int incY){
  int i;
  for(i = 0; i < N; i++, X += incX * binned_zbnum(fold), Y += incY * binned_zbnum(fold)){
    binned_zmzmadd(fold, X, 1, X + 2 * fold, 1, Y, 1, Y + 2 * fold, 1);
  }
}
