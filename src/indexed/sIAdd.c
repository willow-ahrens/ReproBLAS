/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"
#include "../types.h"

// ADDING TWO INDEXED FP
// X += Y
void sIAdd1(int n, float* x, float* xc, int incx, float* y, float* yc, int incy) {
    int i;
    int shift;
    float *repX = y;
    int increpX = incy;
    float *carX = yc;
    int inccarX = incy;
    float *repY = x;
    int increpY = incx;
    float *carY = xc;
    int inccarY = incx;
    int fold = n;

	if (repX[0] == 0.0)
		return;

	if (repY[0] == 0.0) {
		for (i = 0; i < fold; i++) {
			repY[i*increpY] = repX[i*increpX];
			carY[i*inccarY] = carX[i*inccarX];
		}
		return;
	}

    shift = siindex(repY) - siindex(repX);
    if(shift > 0){
      //shift Y upwards and add X
      for (i = fold - 1; i >= shift, i >= 0; i--) {
        repY[i*increpY] = repX[i*increpX] + (repY[(i - shift)*increpY] - 1.5*ufpf(repY[(i - shift)*increpY]));
        carY[i*inccarY] = carX[i*inccarX] + carY[(i - shift)*inccarY];
      }
      for (i = 0; i < shift; i++) {
        repY[i*increpY] = repX[i*increpX];
        carY[i*inccarY] = carX[i*inccarX];
      }
    }else{
      //shift X upwards and add X
	  for (i = 0 - shift; i < fold; i++) {
		repY[i*increpY] += repX[(i + shift)*increpX] - 1.5*ufpf(repX[(i + shift)*increpX]);
		carY[i*inccarY] += carX[(i + shift)*inccarX];
	  }
    }
}
void sIAdd(I_float* X, I_float Y) {
	sIAdd1(DEFAULT_FOLD, (X)->m, (X)->c, 1, (Y).m, (Y).c, 1);
    sirenorm(X, DEFAULT_FOLD);
}

void cIAdd1(int K, float complex* x, float* xc, int incx,
	float complex* y, float* yc, int incy) {
	sIAdd1(K, (float*)x    , xc    , 2 * incx, (float*)y    , yc    , 2 * incy);
	sIAdd1(K, ((float*)x) + 1, xc + 1, 2 * incx, ((float*)y) + 1, yc + 1, 2 * incy);
}

void cIAdd(I_float_Complex* X, I_float_Complex Y) {
	cIAdd1(DEFAULT_FOLD, (float complex*)((X)->m), ((X)->c), 1,
		(float complex*)((Y).m), ((Y).c), 1);	
    cirenorm(X, DEFAULT_FOLD);
}

// no update
void sIAddf1(int fold, float* x, int inc, float y) {
	float M;
	int i;
	int_float iM;
	for (i = 0; i < fold; i++, x += inc) {
		M = x[0];
		iM.f = y;
		iM.i |= 1;
		iM.f += M;
		x[0] = iM.f;
		M    -= iM.f;
		y    += M;
	}
}

void sIAddf(I_float* X, float Y) {
    //printf("X'[0] %g Y %g\n", X->m[0], fabs(Y));
	sisupdate(fabs(Y), X, DEFAULT_FOLD);
    //printf("X[0] %g\n", X->m[0]);
	sIAddf1(DEFAULT_FOLD, (X)->m, 1, Y);
    sirenorm(X, DEFAULT_FOLD);
}

void cIAddc1(int fold, float complex* x, int inc, float complex Y) {
	float* yptr = (float*) &Y;
	sIAddf1(fold, (float*)x, 2 * inc,yptr[0]);
	sIAddf1(fold, ((float*)x)+1, 2*inc,yptr[1]);
}

void cIAddc(I_float_Complex* X, float complex Y) {
    CSET_(Y, fabs(CREAL_(Y)), fabs(CIMAG_(Y)));
    cicupdate(&Y, X, DEFAULT_FOLD);
    //cisupdate(fabs(Y), X, DEFAULT_FOLD);
	cIAddc1(DEFAULT_FOLD, (float complex*)X, 1, Y);
    cirenorm(X, DEFAULT_FOLD);
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

