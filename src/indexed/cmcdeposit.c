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
 * @param priY Y's primary vector
 * @param incpriY stride within Y's primary vector (use every incpriY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   10 Jun 2015
 */
void cmcdeposit(const int fold, const void *X, float *priY, const int incpriY){
  float MR, MI;
  int_float qR, qI;
  int i;
  float xR = ((float*)X)[0];
  float xI = ((float*)X)[1];

  if (ISNANINFF(xR) || ISNANINFF(priY[0])){
    priY[0] += xR;
    smsdeposit(fold, xI, priY + 1, 2 * incpriY);
    return;
  }
  if (ISNANINFF(xI) || ISNANINFF(priY[1])){
    priY[1] += xI;
    smsdeposit(fold, xR, priY, 2 * incpriY);
    return;
  }

  if(smindex0(priY) || smindex0(priY + 1)){
    smsdeposit(fold, xR, priY, 2 * incpriY);
    smsdeposit(fold, xI, priY + 1, 2 * incpriY);
    return;
  }

  for (i = 0; i < fold - 1; i++) {
    MR = priY[i * 2 * incpriY];
    MI = priY[i * 2 * incpriY + 1];
    qR.f = xR;
    qI.f = xI;
    qR.i |= 1;
    qI.i |= 1;
    qR.f += MR;
    qI.f += MI;
    priY[i * 2 * incpriY] = qR.f;
    priY[i * 2 * incpriY + 1] = qI.f;
    MR -= qR.f;
    MI -= qI.f;
    xR += MR;
    xI += MI;
  }
  qR.f = xR;
  qI.f = xI;
  qR.i |= 1;
  qI.i |= 1;
  priY[i * 2 * incpriY] += qR.f;
  priY[i * 2 * incpriY + 1] += qI.f;
}
