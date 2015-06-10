#include <math.h>

#include <indexed.h>

#include "../common/common.h"

/**
 * @internal
 * @brief  Add complex single precision to suitably indexed manually specified indexed complex single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called #cmcupdate() on Y with the maximum absolute value of all future elements you wish to deposit in Y, you can call #cmcdeposit() to deposit a maximum of #SIENDURANCE elements into Y before renormalizing Y with #cmrenorm(). After any number of successive calls of #cmcdeposit() on Y, you must renormalize Y with #cmrenorm() before using any other function on Y.
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
void cmcdeposit(const int fold, const void *X, float *manY, const int incmanY){
  float MR, MI;
  int_float qR, qI;
  int i;
  float xR = ((float*)X)[0];
  float xI = ((float*)X)[1];

  if (isinf(xR) || isnan(xR) || isinf(manY[0]) || isnan(manY[0])) {
    manY[0] += xR;
    smsdeposit(fold, xI, manY + 1, 2 * incmanY);
    return;
  }
  if (isinf(xI) || isnan(xI) || isinf(manY[1]) || isnan(manY[1])) {
    manY[1] += xI;
    smsdeposit(fold, xR, manY, 2 * incmanY);
    return;
  }

  if(smindex0(manY) || smindex0(manY + 1)){
    smsdeposit(fold, xR, manY, 2 * incmanY);
    smsdeposit(fold, xI, manY + 1, 2 * incmanY);
    return;
  }

  for (i = 0; i < fold - 1; i++) {
    MR = manY[i * 2 * incmanY];
    MI = manY[i * 2 * incmanY + 1];
    qR.f = xR;
    qI.f = xI;
    qR.i |= 1;
    qI.i |= 1;
    qR.f += MR;
    qI.f += MI;
    manY[i * 2 * incmanY] = qR.f;
    manY[i * 2 * incmanY + 1] = qI.f;
    MR -= qR.f;
    MI -= qI.f;
    xR += MR;
    xI += MI;
  }
  qR.f = xR;
  qI.f = xI;
  qR.i |= 1;
  qI.i |= 1;
  manY[i * 2 * incmanY] += qR.f;
  manY[i * 2 * incmanY + 1] += qI.f;
}
