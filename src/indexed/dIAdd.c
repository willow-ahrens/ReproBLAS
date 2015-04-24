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
// X += Y
void dIAdd1(int K, double* x, double* xc, int incx, double* y, double* yc, int incy) {
	int i;
    int shift;

    double *repX = y;
    int increpX = incy;
    double *carX = yc;
    int inccarX = incy;
    double *repY = x;
    int increpY = incx;
    double *carY = xc;
    int inccarY = incx;
    int fold = K;

	if (repX[0] == 0.0)
		return;

	if (repY[0] == 0.0) {
		for (i = 0; i < fold; i++) {
			repY[i*increpY] = repX[i*increpX];
			carY[i*inccarY] = carX[i*inccarX];
		}
		return;
	}

    shift = diindex(repY) - diindex(repX);
    if(shift > 0){
      //shift Y upwards and add X
      for (i = fold - 1; i >= shift, i >= 0; i--) {
        repY[i*increpY] = repX[i*increpX] + (repY[(i - shift)*increpY] - 1.5*ufp(repY[(i - shift)*increpY]));
        carY[i*inccarY] = carX[i*inccarX] + carY[(i - shift)*inccarY];
      }
      for (i = 0; i < shift; i++) {
        repY[i*increpY] = repX[i*increpX];
        carY[i*inccarY] = carX[i*inccarX];
      }
    }else{
      //shift X upwards and add X
	  for (i = 0 - shift; i < fold; i++) {
		repY[i*increpY] += repX[(i + shift)*increpX] - 1.5*ufp(repX[(i + shift)*increpX]);
		carY[i*inccarY] += carX[(i + shift)*inccarX];
	  }
    }
}

// X += Y
void dIAdd(I_double* X, I_double Y) {
	dIAdd1(DEFAULT_FOLD, X->m, X->c, 1, (Y).m, (Y).c,1);
    direnorm(X, DEFAULT_FOLD);
}

void zIAdd1(int K,
	double complex* x, double complex* xc, int incx,
	double complex* y, double complex* yc, int incy) {
	dIAdd1(K, (double*)x, (double*)xc, 2*incx, (double*)y, (double*)yc, 2*incy);
	dIAdd1(K, ((double*)x) + 1, ((double*)xc) + 1, 2*incx, ((double*)y)+1, ((double*)yc)+1, 2*incy);
}

void zIAdd(I_double_Complex* X, I_double_Complex Y) {
	zIAdd1(DEFAULT_FOLD,(double complex*)(X)->m,(double complex*)(X)->c,1,
		(double complex*)(Y).m,(double complex*)(Y).c,1);
    zirenorm(X, DEFAULT_FOLD);
}

void dIAddd1(int fold, double* x, int inc, double y) {
	double M;
	long_double lM;
	int i;
	for (i = 0; i < fold - 1; i++, x += inc) {
		M = x[0];
		lM.d = y;
		lM.l |= 1;
		lM.d += M;
		x[0] = lM.d;
		M -= lM.d;
		y += M;
	}
	lM.d = y;
	lM.l |= 1;
	x[0] += lM.d;
}

void dIAddd(I_double* X, double Y) {
	didupdate(fabs(Y), X, DEFAULT_FOLD);
	dIAddd1(DEFAULT_FOLD, X->m, 1, Y);
    direnorm(X, DEFAULT_FOLD);
}

void zIAddz1(int fold, double complex* x, int inc, double complex y) {
	double MR, MI;
	long_double lMR, lMI;
	double* xptr = (double*) x;
	double* yptr = (double*) &y;

	int i;
	double yR = yptr[0];
	double yI = yptr[1];

	inc *= 2;
	for (i = 0; i < fold; i++, xptr += inc) {
		MR = xptr[0];
		MI = xptr[1];

		lMR.d = yR;
		lMI.d = yI;

		lMR.l |= 1;
		lMI.l |= 1;

		lMR.d += MR;
		lMI.d += MI;

		xptr[0] = lMR.d;
		xptr[1] = lMI.d;

		MR -= lMR.d;
		MI -= lMI.d;

		yR += MR;
		yI += MI;
	}
}

//TODO use a zizupdate
void zIAddz(I_double_Complex* X, double complex Y) {
	zidupdate(fabs(Y), X, DEFAULT_FOLD);
	zIAddz1(DEFAULT_FOLD, (double complex*)X->m, 1, Y);
    zirenorm(X, DEFAULT_FOLD);
}

void dINeg1(int fold, double* x, double* c, int inc) {
	double M, X;
	int i;
	for (i = 0; i < fold; i++, x += inc, c += inc) {
		X = x[0];
		M = ufp(X);
		x[0] = (3 * M) - X;
		c[0] = -c[0];
	}
}

void zINeg1(int fold, double complex* x, double complex* c, int inc) {
	double MR, MI, BR, BI;
	double* xptr = (double*) x;
	double* cptr = (double*) c;

	int i;

	inc *= 2;
	for (i = 0; i < fold; i++, xptr += inc, cptr += inc) {
		BR = xptr[0];
		BI = xptr[1];

		MR = ufp(BR);
		MI = ufp(BI);

		xptr[0] = (3 * MR) - BR;
		xptr[1] = (3 * MI) - BI;
		cptr[0] = -cptr[0];
		cptr[1] = -cptr[1];
	}
}

