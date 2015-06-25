#include <indexed.h>

/**
 * @brief  Add indexed complex double precision vectors (Y += X)
 *
 * Performs the operation Y += X
 *
 * @param fold the fold of the indexed types
 * @param N vector length
 * @param X indexed vector X
 * @param incX X vector stride (use every incX'th element)
 * @param Y indexed vector Y
 * @param incY Y vector stride (use every incY'th element)
 *
 * @author Peter Ahrens
 * @date   25 Jun 2015
 */
void ziziaddv(const int fold, const int N, const double_complex_indexed *X, const int incX, double_complex_indexed *Y, const int incY){
  int i;
  for(i = 0; i < N; i++, X += incX * zinum(fold), Y += incY * zinum(fold)){
    zmzmadd(fold, X, 1, X + 2 * fold, 1, Y, 1, Y + 2 * fold, 1);
  }
}
