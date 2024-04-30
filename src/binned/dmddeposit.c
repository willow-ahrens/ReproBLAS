#include <math.h>

#include <binned.h>

#include "../common/common.h"

/**
 * @internal
 * @brief  Add double precision to suitably binned manually specified binned double precision (Y += X)
 *
 * Performs the operation Y += X on an binned type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called #binned_dmdupdate() on Y with the maximum absolute value of all future elements you wish to deposit in Y, you can call #binned_dmddeposit() to deposit a maximum of #binned_DBENDURANCE elements into Y before renormalizing Y with #binned_dmrenorm(). After any number of successive calls of #binned_dmddeposit() on Y, you must renormalize Y with #binned_dmrenorm() before using any other function on Y.
 *
 * @param fold the fold of the binned types
 * @param X scalar X
 * @param priY Y's primary vector
 * @param incpriY stride within Y's primary vector (use every incpriY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Willow Ahrens
 * @date   10 Jun 2015
 */
void binned_dmddeposit(const int fold, const double X, double *priY, const int incpriY){
  double M;
  long_double q;
  int i;
  double x = X;

  if(ISNANINF(x) || ISNANINF(priY[0])){
    priY[0] += x;
    return;
  }

  if(binned_dmindex0(priY)){
    M = priY[0];
    q.d = x * binned_DMCOMPRESSION;
    q.l |= 1;
    q.d += M;
    priY[0] = q.d;
    M -= q.d;
    M *= binned_DMEXPANSION * 0.5;
    x += M;
    x += M;
    for (i = 1; i < fold - 1; i++) {
      M = priY[i * incpriY];
      q.d = x;
      q.l |= 1;
      q.d += M;
      priY[i * incpriY] = q.d;
      M -= q.d;
      x += M;
    }
    q.d = x;
    q.l |= 1;
    priY[i * incpriY] += q.d;
  }else{
    for (i = 0; i < fold - 1; i++) {
      M = priY[i * incpriY];
      q.d = x;
      q.l |= 1;
      q.d += M;
      priY[i * incpriY] = q.d;
      M -= q.d;
      x += M;
    }
    q.d = x;
    q.l |= 1;
    priY[i * incpriY] += q.d;
  }
}
