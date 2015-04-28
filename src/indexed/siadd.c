/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

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

void sisiadd(const int fold, float_indexed *X, float_indexed *Y){
  smsmadd(fold, X, 1, X + fold, 1, Y, 1, Y + fold, 1);
}

void cmcmadd(const int fold, float *repX, int increpX, float *carX, int inccarX, float* repY, int increpY, float* carY, int inccarY) {
  smsmadd(fold, repX, 2 * increpX, carX, 2 * inccarX, repY, 2 * increpY, carY, 2 * inccarY);
  smsmadd(fold, repX + 1, 2 * increpX, carX + 1, 2 * inccarX, repY + 1, 2 * increpY, carY + 1, 2 * inccarY);
}

void ciciadd(const int fold, float_complex_indexed *X, float_complex_indexed *Y){
  cmcmadd(fold, X, 1, X + 2 * fold, 1, Y, 1, Y + 2 * fold, 1);
}

void smsdeposit(const int fold, float X, float *repY, int increpY){
  float M;
  int_float q;
  int i;
  for (i = 0; i < fold - 1; i++) {
    M = repY[i * increpY];
    q.f = X;
    q.i |= 1;
    q.f += M;
    repY[i * increpY] = q.f;
    M -= q.f;
    X += M;
  }
  q.f = X;
  q.i |= 1;
  repY[i * increpY] += q.f;
}

void sisdeposit(const int fold, float X, float_indexed *Y){
  smsdeposit(fold, X, Y, 1);
}

void smsadd(const int fold, float X, float *repY, int increpY, float *carY, int inccarY){
  smsupdate(fold, fabsf(X), repY, increpY, carY, inccarY);
  smsdeposit(fold, X, repY, increpY);
  smrenorm(fold, repY, increpY, carY, inccarY);
}

void sisadd(const int fold, float X, float_indexed *Y){
  smsadd(fold, X, Y, 1, Y + fold, 1);
}

void cmcdeposit(const int fold, void *X, float *repY, int increpY){
  float MR, MI;
  int_float qR, qI;
  int i;
  float xR = ((float*)X)[0];
  float xI = ((float*)X)[1];

  increpY *= 2;

  for (i = 0; i < fold - 1; i++) {
    MR = repY[i * increpY];
    MI = repY[i * increpY + 1];
    qR.f = xR;
    qI.f = xI;
    qR.i |= 1;
    qI.i |= 1;
    qR.f += MR;
    qI.f += MI;
    repY[i * increpY] = qR.f;
    repY[i * increpY + 1] = qI.f;
    MR -= qR.f;
    MI -= qI.f;
    xR += MR;
    xI += MI;
  }
  MR = repY[i * increpY];
  MI = repY[i * increpY + 1];
  qR.f = xR;
  qI.f = xI;
  qR.i |= 1;
  qI.i |= 1;
  qR.f += MR;
  qI.f += MI;
  repY[i * increpY] = qR.f;
  repY[i * increpY + 1] = qI.f;
}

void cicdeposit(const int fold, void *X, float_complex_indexed *Y){
  cmcdeposit(fold, X, Y, 1);
}

void cmcadd(const int fold, void *X, float *repY, int increpY, float *carY, int inccarY){
  float aX[2];
  aX[0] = fabsf(((float*)X)[0]);
  aX[1] = fabsf(((float*)X)[1]);
  cmcupdate(fold, aX, repY, increpY, carY, inccarY);
  cmcdeposit(fold, X, repY, increpY);
  cmrenorm(fold, repY, increpY, carY, inccarY);
}

void cicadd(const int fold, void *X, float_complex_indexed *Y){
  cmcadd(fold, X, Y, 1, Y + 2 * fold, 1);
}
