/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <math.h>
#include <float.h>

#include "indexed.h"
#include "../Common/Common.h"

/**
 * @internal
 * @brief  Add manually specified indexed single precision (Y += X)
 *
 * Performs the operation Y += X
 *
 * @param fold the fold of the indexed types
 * @param manX X's mantissa vector
 * @param incmanX stride within X's mantissa vector (use every incmanX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void smsmadd(const int fold, const float *manX, const int incmanX, const float *carX, const int inccarX, float* manY, const int incmanY, float* carY, const int inccarY) {
  int i;
  int shift;

  if (manX[0] == 0.0)
    return;

  if (manY[0] == 0.0) {
    for (i = 0; i < fold; i++) {
      manY[i*incmanY] = manX[i*incmanX];
      carY[i*inccarY] = carX[i*inccarX];
    }
    return;
  }

  if (isinf(manX[0]) || isnan(manX[0]) || isinf(manY[0]) || isnan(manY[0])) {
    manY[0] += manX[0];
    return;
  }

  shift = smindex(manY) - smindex(manX);
  if(shift > 0){
    //shift Y upwards and add X to Y
    for (i = fold - 1; i >= shift; i--) {
      manY[i*incmanY] = manX[i*incmanX] + (manY[(i - shift)*incmanY] - 1.5*ufpf(manY[(i - shift)*incmanY]));
      carY[i*inccarY] = carX[i*inccarX] + carY[(i - shift)*inccarY];
    }
    for (i = 0; i < shift && i < fold; i++) {
      manY[i*incmanY] = manX[i*incmanX];
      carY[i*inccarY] = carX[i*inccarX];
    }
  }else{
    //shift X upwards and add X to Y
    for (i = 0 - shift; i < fold; i++) {
      manY[i*incmanY] += manX[(i + shift)*incmanX] - 1.5*ufpf(manX[(i + shift)*incmanX]);
      carY[i*inccarY] += carX[(i + shift)*inccarX];
    }
  }

  smrenorm(fold, manY, incmanY, carY, inccarY);
}

/**
 * @brief  Add indexed single precision (Y += X)
 *
 * Performs the operation Y += X
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void sisiadd(const int fold, const float_indexed *X, float_indexed *Y){
  smsmadd(fold, X, 1, X + fold, 1, Y, 1, Y + fold, 1);
}

/**
 * @internal
 * @brief  Add manually specified indexed complex single precision (Y += X)
 *
 * Performs the operation Y += X
 *
 * @param fold the fold of the indexed types
 * @param manX X's mantissa vector
 * @param incmanX stride within X's mantissa vector (use every incmanX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cmcmadd(const int fold, const float *manX, const int incmanX, const float *carX, const int inccarX, float* manY, const int incmanY, float* carY, const int inccarY) {
  smsmadd(fold, manX, 2 * incmanX, carX, 2 * inccarX, manY, 2 * incmanY, carY, 2 * inccarY);
  smsmadd(fold, manX + 1, 2 * incmanX, carX + 1, 2 * inccarX, manY + 1, 2 * incmanY, carY + 1, 2 * inccarY);
}

/**
 * @brief  Add indexed complex single precision (Y += X)
 *
 * Performs the operation Y += X
 *
 * @param fold the fold of the indexed types
 * @param X indexed scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void ciciadd(const int fold, const float_complex_indexed *X, float_complex_indexed *Y){
  cmcmadd(fold, X, 1, X + 2 * fold, 1, Y, 1, Y + 2 * fold, 1);
}

/**
 * @internal
 * @brief  Add single precision to suitably indexed manually specified indexed single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called smsupdate() on Y with the maximum absolute value of any elements you wish to deposit in Y, you can call this method to deposit a maximum of sicapacity() elements into Y. After calling smsdeposit() on an indexed type, you must renormalize the indexed type with smrenorm().
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
void smsdeposit(const int fold, const float X, float *manY, const int incmanY){
  float M;
  int_float q;
  int i;
  float x = X * smcompression();

  if (isinf(x) || isnan(x) || isinf(manY[0]) || isnan(manY[0])) {
    manY[0] += x;
    return;
  }

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

/**
 * @brief  Add single precision to suitably indexed indexed single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called sisupdate() on Y with the maximum absolute value of any elements you wish to deposit in Y, you can call this method to deposit a maximum of sicapacity() elements into Y. After calling sisdeposit() on an indexed type, you must renormalize the indexed type with sirenorm().
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void sisdeposit(const int fold, const float X, float_indexed *Y){
  smsdeposit(fold, X, Y, 1);
}

/**
 * @internal
 * @brief  Add complex single precision to suitably indexed manually specified indexed complex single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called cmcupdate() on Y with the maximum absolute value of any elements you wish to deposit in Y, you can call this method to deposit a maximum of sicapacity() elements into Y. After calling cmcdeposit() on an indexed type, you must renormalize the indexed type with cmrenorm().
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
void cmcdeposit(const int fold, const void *X, float *manY, const int incmanY){
  float MR, MI;
  int_float qR, qI;
  int i;
  float xR = ((float*)X)[0] * smcompression();
  float xI = ((float*)X)[1] * smcompression();

  if (isinf(xR) || isnan(xR) || isinf(manY[0]) || isnan(manY[0])) {
    manY[0] += xR;
    smsdeposit(fold, xI, manY + incmanY, 2 * incmanY);
    return;
  }
  if (isinf(xI) || isnan(xI) || isinf(manY[1]) || isnan(manY[1])) {
    manY[1] += xI;
    smsdeposit(fold, xR, manY, 2 * incmanY);
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

/**
 * @brief  Add complex single precision to suitably indexed manually specified indexed complex single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called cicupdate() on Y with the maximum absolute value of any elements you wish to deposit in Y, you can call this method to deposit a maximum of sicapacity() elements into Y. After calling cicdeposit() on an indexed type, you must renormalize the indexed type with cirenorm().
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
void cicdeposit(const int fold, const void *X, float_complex_indexed *Y){
  cmcdeposit(fold, X, Y, 1);
}

/**
 * @internal
 * @brief  Add single precision to manually specified indexed single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void smsadd(const int fold, const float X, float *manY, const int incmanY, float *carY, const int inccarY){
  smsupdate(fold, fabsf(X), manY, incmanY, carY, inccarY);
  smsdeposit(fold, X, manY, incmanY);
  smrenorm(fold, manY, incmanY, carY, inccarY);
}

/**
 * @brief  Add single precision to indexed single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void sisadd(const int fold, const float X, float_indexed *Y){
  smsadd(fold, X, Y, 1, Y + fold, 1);
}

/**
 * @internal
 * @brief  Add complex single precision to manually specified indexed complex single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param manY Y's mantissa vector
 * @param incmanY stride within Y's mantissa vector (use every incmanY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cmcadd(const int fold, const void *X, float *manY, const int incmanY, float *carY, const int inccarY){
  float aX[2];
  aX[0] = fabsf(((float*)X)[0]);
  aX[1] = fabsf(((float*)X)[1]);
  cmcupdate(fold, aX, manY, incmanY, carY, inccarY);
  cmcdeposit(fold, X, manY, incmanY);
  cmrenorm(fold, manY, incmanY, carY, inccarY);
}

/**
 * @brief  Add complex single precision to indexed complex single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param Y indexed scalar Y
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cicadd(const int fold, const void *X, float_complex_indexed *Y){
  cmcadd(fold, X, Y, 1, Y + 2 * fold, 1);
}
