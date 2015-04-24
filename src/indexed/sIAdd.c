/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

void smsmadd(float *repX, int increpX, float *carX, int inccarX, float* repY, int increpY, float* carY, int inccarY, int fold) {
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

  smrenorm(repY, increpY, carY, inccarY, fold);
}

void sisiadd(float_indexed *X, float_indexed *Y, int fold){
  smsmadd(X, 1, X + fold, 1, Y, 1, Y + fold, 1, fold);
}

void cmcmadd(float *repX, int increpX, float *carX, int inccarX, float* repY, int increpY, float* carY, int inccarY, int fold) {
  smsmadd(repX, 2 * increpX, carX, 2 * inccarX, repY, 2 * increpY, carY, 2 * inccarY, fold);
  smsmadd(repX + 1, 2 * increpX, carX + 1, 2 * inccarX, repY + 1, 2 * increpY, carY + 1, 2 * inccarY, fold);
}

void ciciadd(float_complex_indexed *X, float_complex_indexed *Y, int fold){
  cmcmadd(X, 1, X + 2 * fold, 1, Y, 1, Y + 2 * fold, 1, fold);
}

void smsdeposit(float X, float *repY, int increpY, int fold){
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

void sisdeposit(float X, float_indexed *Y, int fold){
  smsdeposit(X, Y, 1, fold);
}

void smsadd(float X, float *repY, int increpY, float *carY, int inccarY, int fold){
  smsupdate(fabsf(X), repY, increpY, carY, inccarY, fold);
  smsdeposit(X, repY, increpY, fold);
  smrenorm(repY, increpY, carY, inccarY, fold);
}

void sisadd(float X, float_indexed *Y, int fold){
  smsadd(X, Y, 1, Y + fold, 1, fold);
}

void cmcdeposit(void *X, float *repY, int increpY, fold){
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

void cicdeposit(void *X, float_complex_indexed *Y, fold){
  cmcdeposit(X, Y, 1, fold);
}

void cmcadd(void *X, float *repY, int increpY, float *carY, int inccarY, int fold){
  float aX[2];
  aX[0] = fabsf(((float*)X)[0]);
  aX[1] = fabsf(((float*)X)[1]);
  cmcupdate(aX, repY, increpY, carY, inccarY, fold);
  cmcdeposit(X, repY, increpY, fold);
  cmrenorm(repY, increpY, carY, inccarY, fold);
}

void cicadd(void *X, float_complex_indexed *Y, int fold){
  cmcadd(X, Y, 1, Y + 2 * fold, 1, fold);
}

void sINeg1(int fold, float* x, float* c, int inc) {
	float M, X;
	int i;
	for (i = 0; i < fold; i++, x += inc, c += inc) {
		X = x[0];
		M = ufpf(X);
		x[0] = (3 * M) - X;
		c[0] = -c[0];
	}
}

void cINeg1(int fold, float complex* x, float* c, int inc) {
	float MR, MI, BR, BI;
	float* xptr = (float*) x;
	int i;

	inc *= 2;
	for (i = 0; i < fold; i++, xptr += inc, c += inc) {
		BR = xptr[0];
		BI = xptr[1];

		MR = ufpf(BR);
		MI = ufpf(BI);

		xptr[0] = (3 * MR) - BR;
		xptr[1] = (3 * MI) - BI;
		c[0] = -c[0];
		c[1] = -c[1];
	}
}

