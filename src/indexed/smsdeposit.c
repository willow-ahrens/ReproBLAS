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
 * @param priY Y's primary vector
 * @param incpriY stride within Y's primary vector (use every incpriY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   10 Jun 2015
 */
void smsdeposit(const int fold, const float X, float *priY, const int incpriY){
  float M;
  int_float q;
  int i;
  float x = X;

  if (ISNANINFF(x) || ISNANINFF(priY[0])) {
    priY[0] += x;
    return;
  }

  if(smindex0(priY)){
    M = priY[0];
    q.f = x * SMCOMPRESSION;
    q.i |= 1;
    q.f += M;
    priY[0] = q.f;
    if (fold > 1) {
      M -= q.f;
      x += M * SMEXPANSION;
      for (i = 1; i < fold - 1; i++) {
        M = priY[i * incpriY];
        q.f = x;
        q.i |= 1;
        q.f += M;
        priY[i * incpriY] = q.f;
        M -= q.f;
        x += M;
      }
      q.f = x;
      q.i |= 1;
      priY[i * incpriY] += q.f;
    }
  }else{
    for (i = 0; i < fold - 1; i++) {
      M = priY[i * incpriY];
      q.f = x;
      q.i |= 1;
      q.f += M;
      priY[i * incpriY] = q.f;
      M -= q.f;
      x += M;
    }
    q.f = x;
    q.i |= 1;
    priY[i * incpriY] += q.f;
  }
}
