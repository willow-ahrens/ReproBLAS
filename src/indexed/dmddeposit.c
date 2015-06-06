#include <math.h>

#include <indexed.h>

#include "../common/common.h"

/**
 * @internal
 * @brief  Add double precision to suitably indexed manually specified indexed double precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called dmdupdate() on Y with the maximum absolute value of any elements you wish to deposit in Y, you can call this method to deposit a maximum of DIENDURANCE elements into Y. After calling dmddeposit() on an indexed type, you must renormalize the indexed type with dmrenorm().
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void dmddeposit(const int fold, const double X, double *manY, const int incmanY){
  double M;
  long_double q;
  int i;
  double x = X;

  if (isinf(x) || isnan(x) || isinf(manY[0]) || isnan(manY[0])) {
    manY[0] += x;
    return;
  }

  if(dmindex0(manY)){
    M = manY[0];
    q.d = x * DMCOMPRESSION;
    q.l |= 1;
    q.d += M;
    manY[0] = q.d;
    if (fold > 1) {
      M -= q.d;
      x += M * DMEXPANSION;
      for (i = 1; i < fold - 1; i++) {
        M = manY[i * incmanY];
        q.d = x;
        q.l |= 1;
        q.d += M;
        manY[i * incmanY] = q.d;
        M -= q.d;
        x += M;
      }
      q.d = x;
      q.l |= 1;
      manY[i * incmanY] += q.d;
    }
  }else{
    for (i = 0; i < fold - 1; i++) {
      M = manY[i * incmanY];
      q.d = x;
      q.l |= 1;
      q.d += M;
      manY[i * incmanY] = q.d;
      M -= q.d;
      x += M;
    }
    q.d = x;
    q.l |= 1;
    manY[i * incmanY] += q.d;
  }
}
