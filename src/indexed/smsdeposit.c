#include <math.h>

#include <indexed.h>

#include "../common/common.h"

/**
 * @internal
 * @brief  Add single precision to suitably indexed manually specified indexed single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called #smsupdate() on Y with the maximum absolute value of all future elements you wish to deposit in Y, you can call #smsdeposit() to deposit a maximum of #SIENDURANCE elements into Y before renormalizing Y with #smrenorm(). After any number of successive calls of #smsdeposit() on Y, you must renormalize Y with #smrenorm() before using any other function on Y.
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   10 Jun 2015
 */
void smsdeposit(const int fold, const float X, float *manY, const int incmanY){
  float M;
  int_float q;
  int i;
  float x = X;

  if (ISNANINFF(x) || ISNANINFF(manY[0])) {
    manY[0] += x;
    return;
  }

  if(smindex0(manY)){
    M = manY[0];
    q.f = x * SMCOMPRESSION;
    q.i |= 1;
    q.f += M;
    manY[0] = q.f;
    if (fold > 1) {
      M -= q.f;
      x += M * SMEXPANSION;
      for (i = 1; i < fold - 1; i++) {
        M = manY[i * incmanY];
        q.f = x;
        q.i |= 1;
        q.f += M;
        manY[i * incmanY] = q.f;
        M -= q.f;
        x += M;
      }
      q.f = x;
      q.i |= 1;
      manY[i * incmanY] += q.f;
    }
  }else{
    for (i = 0; i < fold - 1; i++) {
      M = manY[i * incmanY];
      q.f = x;
      q.i |= 1;
      q.f += M;
      manY[i * incmanY] = q.f;
      M -= q.f;
      x += M;
    }
    q.f = x;
    q.i |= 1;
    manY[i * incmanY] += q.f;
  }
}
