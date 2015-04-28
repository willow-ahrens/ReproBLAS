/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

/**
 * @brief  Add manually specified indexed single precision (Y += X)
 *
 * Performs the operation Y += X
 *
 * @param fold the fold of the indexed types
 * @param repX X's rep vector
 * @param increpX stride within X's rep vector (use every increpX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param repY Y's rep vector
 * @param increpY stride within Y's rep vector (use every increpY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void smsmadd(const int fold, float *repX, int increpX, float *carX, int inccarX, float* repY, int increpY, float* carY, int inccarY) {
  int i;
  int shift;

  if (repX[0] == 0.0)
    return;

  if (repY[0] == 0.0) {
    for (i = 0; i < fold; i++) {
      repY[i*increpY] = repX[i*increpX];
      carY[i*inccarY] = carX[i*inccarX];
    }
    return;
  }

  shift = smindex(repY) - smindex(repX);
  if(shift > 0){
    //shift Y upwards and add X to Y
    for (i = fold - 1; i >= shift; i--) {
      repY[i*increpY] = repX[i*increpX] + (repY[(i - shift)*increpY] - 1.5*ufpf(repY[(i - shift)*increpY]));
      carY[i*inccarY] = carX[i*inccarX] + carY[(i - shift)*inccarY];
    }
    for (i = 0; i < shift && i < fold; i++) {
      repY[i*increpY] = repX[i*increpX];
      carY[i*inccarY] = carX[i*inccarX];
    }
  }else{
    //shift X upwards and add X to Y
    for (i = 0 - shift; i < fold; i++) {
      repY[i*increpY] += repX[(i + shift)*increpX] - 1.5*ufpf(repX[(i + shift)*increpX]);
      carY[i*inccarY] += carX[(i + shift)*inccarX];
    }
  }

  smrenorm(fold, repY, increpY, carY, inccarY);
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
void sisiadd(const int fold, float_indexed *X, float_indexed *Y){
  smsmadd(fold, X, 1, X + fold, 1, Y, 1, Y + fold, 1);
}

/**
 * @brief  Add manually specified indexed complex single precision (Y += X)
 *
 * Performs the operation Y += X
 *
 * @param fold the fold of the indexed types
 * @param repX X's rep vector
 * @param increpX stride within X's rep vector (use every increpX'th element)
 * @param carX X's carry vector
 * @param inccarX stride within X's carry vector (use every inccarX'th element)
 * @param repY Y's rep vector
 * @param increpY stride within Y's rep vector (use every increpY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cmcmadd(const int fold, float *repX, int increpX, float *carX, int inccarX, float* repY, int increpY, float* carY, int inccarY) {
  smsmadd(fold, repX, 2 * increpX, carX, 2 * inccarX, repY, 2 * increpY, carY, 2 * inccarY);
  smsmadd(fold, repX + 1, 2 * increpX, carX + 1, 2 * inccarX, repY + 1, 2 * increpY, carY + 1, 2 * inccarY);
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
void ciciadd(const int fold, float_complex_indexed *X, float_complex_indexed *Y){
  cmcmadd(fold, X, 1, X + 2 * fold, 1, Y, 1, Y + 2 * fold, 1);
}

/**
 * @brief  Add single precision to suitably indexed manually specified indexed single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called smsupdate() on Y with the maximum absolute value of any elements you wish to deposit in Y, you can call this method to deposit a maximum of sicapacity() elements into Y. After calling smsdeposit() on an indexed type, you must renormalize the indexed type with smrenorm().
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param repY Y's rep vector
 * @param increpY stride within Y's rep vector (use every increpY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void smsdeposit(const int fold, float X, float *repY, int increpY){
  float M;
  int_float q;
  int i;
  float x = X;
  for (i = 0; i < fold - 1; i++) {
    M = repY[i * increpY];
    q.f = x;
    q.i |= 1;
    q.f += M;
    repY[i * increpY] = q.f;
    M -= q.f;
    x += M;
  }
  q.f = x;
  q.i |= 1;
  repY[i * increpY] += q.f;
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
void sisdeposit(const int fold, float X, float_indexed *Y){
  smsdeposit(fold, X, Y, 1);
}

/**
 * @brief  Add complex single precision to suitably indexed manually specified indexed complex single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y where the index of Y is larger than the index of X
 *
 * @note This routine was provided as a means of allowing the you to optimize your code. After you have called cmcupdate() on Y with the maximum absolute value of any elements you wish to deposit in Y, you can call this method to deposit a maximum of sicapacity() elements into Y. After calling cmcdeposit() on an indexed type, you must renormalize the indexed type with cmrenorm().
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param repY Y's rep vector
 * @param increpY stride within Y's rep vector (use every increpY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cmcdeposit(const int fold, void *X, float *repY, int increpY){
  float MR, MI;
  int_float qR, qI;
  int i;
  float xR = ((float*)X)[0];
  float xI = ((float*)X)[1];

  for (i = 0; i < fold - 1; i++) {
    MR = repY[i * 2 * increpY];
    MI = repY[i * 2 * increpY + 1];
    qR.f = xR;
    qI.f = xI;
    qR.i |= 1;
    qI.i |= 1;
    qR.f += MR;
    qI.f += MI;
    repY[i * 2 * increpY] = qR.f;
    repY[i * 2 * increpY + 1] = qI.f;
    MR -= qR.f;
    MI -= qI.f;
    xR += MR;
    xI += MI;
  }
  MR = repY[i * 2 * increpY];
  MI = repY[i * 2 * increpY + 1];
  qR.f = xR;
  qI.f = xI;
  qR.i |= 1;
  qI.i |= 1;
  qR.f += MR;
  qI.f += MI;
  repY[i * 2 * increpY] = qR.f;
  repY[i * 2 * increpY + 1] = qI.f;
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
 * @param repY Y's rep vector
 * @param increpY stride within Y's rep vector (use every increpY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cicdeposit(const int fold, void *X, float_complex_indexed *Y){
  cmcdeposit(fold, X, Y, 1);
}

/**
 * @brief  Add single precision to manually specified indexed single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param repY Y's rep vector
 * @param increpY stride within Y's rep vector (use every increpY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void smsadd(const int fold, float X, float *repY, int increpY, float *carY, int inccarY){
  smsupdate(fold, fabsf(X), repY, increpY, carY, inccarY);
  smsdeposit(fold, X, repY, increpY);
  smrenorm(fold, repY, increpY, carY, inccarY);
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
void sisadd(const int fold, float X, float_indexed *Y){
  smsadd(fold, X, Y, 1, Y + fold, 1);
}

/**
 * @brief  Add complex single precision to manually specified indexed complex single precision (Y += X)
 *
 * Performs the operation Y += X on an indexed type Y
 *
 * @param fold the fold of the indexed types
 * @param X scalar X
 * @param repY Y's rep vector
 * @param increpY stride within Y's rep vector (use every increpY'th element)
 * @param carY Y's carry vector
 * @param inccarY stride within Y's carry vector (use every inccarY'th element)
 *
 * @author Hong Diep Nguyen
 * @author Peter Ahrens
 * @date   27 Apr 2015
 */
void cmcadd(const int fold, void *X, float *repY, int increpY, float *carY, int inccarY){
  float aX[2];
  aX[0] = fabsf(((float*)X)[0]);
  aX[1] = fabsf(((float*)X)[1]);
  cmcupdate(fold, aX, repY, increpY, carY, inccarY);
  cmcdeposit(fold, X, repY, increpY);
  cmrenorm(fold, repY, increpY, carY, inccarY);
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
void cicadd(const int fold, void *X, float_complex_indexed *Y){
  cmcadd(fold, X, Y, 1, Y + 2 * fold, 1);
}
