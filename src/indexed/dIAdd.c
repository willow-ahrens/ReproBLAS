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
	int i, j;
	int d;

	// Y  == 0
	if (y[0] == 0.0)
		return;

	// X  == 0
	if (x[0] == 0.0) {
		for (i = 0; i < K; i++) {
			x[i*incx] = y[i*incy];
			xc[i*incx] = yc[i*incy];
		}
		return;
	}

	// X AND Y HAVE THE SAME INDEX, JUST ADDING THE CORRESPONDING COMPONENTS
	double MX, MY;
	MX = ufp(x[0]);
	MY = ufp(y[0]);
	if (MX == MY) {
		x[0] += (y[0] - 1.5 * MX);
		for (i = 1; i < K; i++) {
			MX = ufp(x[i*incx]);
			x[i*incx] += (y[i*incy] - 1.5 * MX);
		}
		for (i = 0; i < K; i++)
			xc[i*incx] += yc[i*incy];
		return;
	}

	// INDEX(X) > INDEX(Y): RIGHT-SHIFT Y BEFORE ADDING
	if (MY < MX) {
		d = 0;
		while (d < K - 1 && MY < MX) MX = ufp(x[(++d) * incx]);
		if (MY < MX) return;
		for (i = d, j = 0; i < K; i++, j++) {
			MX = ufp(x[i*incx]);
			x[i*incx]  += (y[j*incy] - 1.5*MX);
			xc[i*incx] += yc[j*incy];
		}
		return;
	}

	// INDEX(X) < INDEX(Y): SHIFT RIGHT X
	d = 0;
	while (d < K - 1 && MX < MY) {
		MY = ufp(y[++d * incy]);
	}
	if (MX < MY) d = K;
	for (i = K - 1; i >= d; i--) {
		MX = ufp(y[i*incy]);
		x[i*incx] = y[i*incy] + (x[(i - d) * incx] - 1.5*MX);
		xc[i*incx] = yc[i*incy] + xc[(i - d)*incx];
	}
	for (i = 0; i < d; i++) {
		x[i*incx] = y[i*incy];
		xc[i*incx] = yc[i*incy];
	}
}

// X += Y
void dIAdd(I_double* X, I_double Y) {
	dIAdd1(DEFAULT_FOLD, X->m, X->c, 1, (Y).m, (Y).c,1);
	dIRenorm1(DEFAULT_FOLD, (X)->m, (X)->c,1);
}

void zIAdd1(int K,
	double complex* x, double complex* xc, int incx,
	double complex* y, double complex* yc, int incy) {
	dIAdd1(K, (double*)x, (double*)xc, 2*incx, (double*)y, (double*)yc, 2*incy);
	dIAdd1(K, (double*)x + 1, (double*)xc + 1, 2*incx, (double*)y+1, (double*)yc+1, 2*incy);
}

void zIAdd(I_double_Complex* X, I_double_Complex Y) {
	zIAdd1(DEFAULT_FOLD,(double complex*)(X)->m,(double complex*)(X)->c,1,
		(double complex*)(Y).m,(double complex*)(Y).c,1);
	zIRenorm1(DEFAULT_FOLD,(double complex*)(X)->m,(double complex*)(X)->c,1);
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
	dIUpdate1(DEFAULT_FOLD, X->m, X->c, 1, fabs(Y));
	dIAddd1(DEFAULT_FOLD, X->m, 1, Y);
	dIRenorm1(DEFAULT_FOLD, X->m, X->c, 1);
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

void zIAddz(I_double_Complex* X, double complex Y) {
	zIUpdate1(DEFAULT_FOLD, 
		(double complex*)X->m, (double complex*)X->c, 1, fabs(Y));
	zIAddz1(DEFAULT_FOLD, (double complex*)X->m, 1, Y);
	zIRenorm1(DEFAULT_FOLD, (double complex*)X->m, (double complex*)X->c, 1);
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

