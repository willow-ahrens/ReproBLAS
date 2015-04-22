/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

// ADDING TWO INDEXED FP
// X += Y
void sIAdd1(int n, float* x, float* xc, int incx, float* y, float* yc, int incy) {
	int i, j;
	int d;

	// Y  == 0
	if (y[0] == 0.0)
		return;

	// X  == 0
	if (x[0] == 0.0) {
		for (i = 0; i < n; i++) {
			x[i*incx]  = y[i*incy];
			xc[i*incx] = yc[i*incy];
		}
		return;
	}

	// X AND Y HAVE THE SAME INDEX, JUST ADDING THE CORRESPONDING COMPONENTS
	float MX, MY;
	MX = ufpf(x[0]);
	MY = ufpf(y[0]);
	if (MX == MY) {
		x[0]  += (y[0] - 1.5 * MX);
		xc[0] += yc[0];
		for (i = 1; i < n; i++) {
			MX = ufpf(x[i*incx]);
			x[i*incx] += (y[i*incy] - 1.5 * MX);
			xc[i*incx] += yc[i*incy];
		}
		return;
	}
	// INDEX(X) > INDEX(Y): RIGHT-SHIFT Y BEFORE ADDING
	if (MY < MX) {
		d = 0;
		while (d < n - 1 && MY < MX) MX = ufpf(x[(++d) * incx]);
		if (MY < MX) return;
		for (i = d, j = 0; i < n; i++, j++) {
			MX = ufpf(x[i*incx]);
			x[i*incx] += (y[j*incy] - 1.5*MX);
			xc[(i )*incx] += yc[j *incy];
		}
		return;
	}

	// INDEX(X) < INDEX(Y): SHIFT RIGHT X
	d = 0;
	while (d < n - 1 && MX < MY) {
		MY = ufpf(y[++d * incy]);
	}
	if (MX < MY) d = n;
	for (i = n - 1; i >= d; i--) {
		MX = ufpf(y[i*incy]);
		x[i*incx] = y[i*incy] + (x[(i - d) * incx] - 1.5*MX);
		xc[(i)*incx] = yc[(i)*incy] + xc[(i- d)*incx];
	}
	for (i = 0; i < d; i++) {
		x[i*incx] = y[i*incy];
		xc[(i)*incx] = yc[(i)*incy];
	}
}
void sIAdd(I_float* X, I_float Y) {
	sIAdd1(DEFAULT_FOLD, (X)->m, (X)->c, 1, (Y).m, (Y).c, 1);
	sIRenorm1(DEFAULT_FOLD,(X)->m,(X)->c,1);
}

void cIAdd1(int K, float complex* x, float* xc, int incx,
	float complex* y, float* yc, int incy) {
	sIAdd1(K, (float*)x    , xc    , 2 * incx, (float*)y    , yc    , 2 * incy);
	sIAdd1(K, (float*)x + 1, xc + 1, 2 * incx, (float*)y + 1, yc + 1, 2 * incy);
}

void cIAdd(I_float_Complex* X, I_float_Complex Y) {
	cIAdd1(DEFAULT_FOLD, (float complex*)(X)->m, (X)->c, 1,
		(float complex*)(Y).m, (Y).c, 1);	
	cIRenorm1(DEFAULT_FOLD, (float complex*)(X)->m, (X)->c, 1);
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
	sisupdate(fabs(Y), X, DEFAULT_FOLD);
	sIAddf1(DEFAULT_FOLD, (X)->m, 1, Y);
	sIRenorm1(DEFAULT_FOLD, (X)->m, (X)->c, 1);
}

void cIAddc1(int fold, float complex* x, int inc, float complex Y) {
	float* yptr = (float*) &Y;
	sIAddf1(fold, (float*)x, 2 * inc,yptr[0]);
	sIAddf1(fold, (float*)x+1, 2*inc,yptr[1]);
}

//TODO use a cicupdate
void cIAddc(I_float_Complex* X, float complex Y) {
	cisupdate(fabs(Y), X, DEFAULT_FOLD);
	cIAddc1(DEFAULT_FOLD, (float complex*)((X)->m), 1, Y);
	cIRenorm1(DEFAULT_FOLD,(float complex*)((X)->m), (X)->c, 1);
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

