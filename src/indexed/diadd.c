/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include "indexed.h"
#include "../Common/Common.h"

// ADDING TWO INDEXED FP
// Y += X
void dmdmadd(double *repX, int increpX, double *carX, int inccarX, double* repY, int increpY, double* carY, int inccarY, int fold) {
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

  shift = dmindex(repY) - dmindex(repX);
  if(shift > 0){
    //shift Y upwards and add X to Y
    for (i = fold - 1; i >= shift; i--) {
      repY[i*increpY] = repX[i*increpX] + (repY[(i - shift)*increpY] - 1.5*ufp(repY[(i - shift)*increpY]));
      carY[i*inccarY] = carX[i*inccarX] + carY[(i - shift)*inccarY];
    }
    for (i = 0; i < shift && i < fold; i++) {
      repY[i*increpY] = repX[i*increpX];
      carY[i*inccarY] = carX[i*inccarX];
    }
  }else{
    //shift X upwards and add X to Y
    for (i = 0 - shift; i < fold; i++) {
      repY[i*increpY] += repX[(i + shift)*increpX] - 1.5*ufp(repX[(i + shift)*increpX]);
      carY[i*inccarY] += carX[(i + shift)*inccarX];
    }
  }

  dmrenorm(repY, increpY, carY, inccarY, fold);
}

void didiadd(double_indexed *X, double_indexed *Y, int fold){
  dmdmadd(X, 1, X + fold, 1, Y, 1, Y + fold, 1, fold);
}

void zmzmadd(double *repX, int increpX, double *carX, int inccarX, double* repY, int increpY, double* carY, int inccarY, int fold) {
  dmdmadd(repX, 2 * increpX, carX, 2 * inccarX, repY, 2 * increpY, carY, 2 * inccarY, fold);
  dmdmadd(repX + 1, 2 * increpX, carX + 1, 2 * inccarX, repY + 1, 2 * increpY, carY + 1, 2 * inccarY, fold);
}

void ziziadd(double_complex_indexed *X, double_complex_indexed *Y, int fold){
  zmzmadd(X, 1, X + 2 * fold, 1, Y, 1, Y + 2 * fold, 1, fold);
}

void dmddeposit(double X, double *repY, int increpY, int fold){
  double M;
  long_double q;
  int i;
  for (i = 0; i < fold - 1; i++) {
    M = repY[i * increpY];
    q.d = X;
    q.l |= 1;
    q.d += M;
    repY[i * increpY] = q.d;
    M -= q.d;
    X += M;
  }
  q.d = X;
  q.l |= 1;
  repY[i * increpY] += q.d;
}

void diddeposit(double X, double_indexed *Y, int fold){
  dmddeposit(X, Y, 1, fold);
}

void dmdadd(double X, double *repY, int increpY, double *carY, int inccarY, int fold){
  dmdupdate(fabs(X), repY, increpY, carY, inccarY, fold);
  dmddeposit(X, repY, increpY, fold);
  dmrenorm(repY, increpY, carY, inccarY, fold);
}

void didadd(double X, double_indexed *Y, int fold){
  dmdadd(X, Y, 1, Y + fold, 1, fold);
}

void zmzdeposit(void *X, double *repY, int increpY, int fold){
  double MR, MI;
  long_double qR, qI;
  int i;
  double xR = ((double*)X)[0];
  double xI = ((double*)X)[1];

  increpY *= 2;

  for (i = 0; i < fold - 1; i++) {
    MR = repY[i * increpY];
    MI = repY[i * increpY + 1];
    qR.d = xR;
    qI.d = xI;
    qR.l |= 1;
    qI.l |= 1;
    qR.d += MR;
    qI.d += MI;
    repY[i * increpY] = qR.d;
    repY[i * increpY + 1] = qI.d;
    MR -= qR.d;
    MI -= qI.d;
    xR += MR;
    xI += MI;
  }
  MR = repY[i * increpY];
  MI = repY[i * increpY + 1];
  qR.d = xR;
  qI.d = xI;
  qR.l |= 1;
  qI.l |= 1;
  qR.d += MR;
  qI.d += MI;
  repY[i * increpY] = qR.d;
  repY[i * increpY + 1] = qI.d;
}

void zizdeposit(void *X, double_complex_indexed *Y, int fold){
  zmzdeposit(X, Y, 1, fold);
}

void zmzadd(void *X, double *repY, int increpY, double *carY, int inccarY, int fold){
  double aX[2];
  aX[0] = fabs(((double*)X)[0]);
  aX[1] = fabs(((double*)X)[1]);
  zmzupdate(aX, repY, increpY, carY, inccarY, fold);
  zmzdeposit(X, repY, increpY, fold);
  zmrenorm(repY, increpY, carY, inccarY, fold);
}

void zizadd(void *X, double_complex_indexed *Y, int fold){
  zmzadd(X, Y, 1, Y + 2 * fold, 1, fold);
}
