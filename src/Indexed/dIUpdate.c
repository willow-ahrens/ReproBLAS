/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "Indexed.h"
#include "../Common/Common.h"

#ifdef __SSE2__
#	include <emmintrin.h>
#endif

#define PREC 53

// max < 2^(-53) * ufp(M) * 2^W

// CHANGE THE INDEX
// IN FACT, IT IS A POSSIBLE SHIFTING LEFT
void dIUpdate1_(int fold, int W, int NB, double step,
	double* x, double* c, int ldx, double y) {
	int i;
	int shift;

	if (y == 0)
		return;

	if (x[0] == 0) {
		dIBoundary(fold, W, fabs(y), x, ldx);
		for (i = fold; i < fold; i++) {
			x[i * ldx] = 0;
			c[i * ldx] = 0;
		}
		return;
	}


	double M = ufp(x[0]);
//	printf("x=%g, M = %g.", x[0], M);

	M = M / NB;

	// NO NEED TO UPDATE
	if (y < M)
		return;

	step = 1.0/step;
	
	int d = 0;
	while (y >= M && d < fold) { d++; M *= step; }

//	printf("Update to new M: %g\n", M);

	if (d == fold) {
		dIBoundary(fold, W, fabs(y), x, ldx);
		for (i = 0; i < fold; i++) {
			c[i*ldx] = 0;
		}
		return;
	}

	for (i = fold-1; i >= d; i--) {
		x[i * ldx] = x[(i - d) * ldx];
		c[i * ldx] = c[(i - d) * ldx];
	}
	dIBoundary(d, W, fabs(y), x, ldx);
	for (i = 0; i < d; i++) c[i * ldx] = 0.0;
}

void dIUpdate1(int fold, int W, double* x, double* c, int ldx, double y) {
	if (y == 0)
		return;

	if (isnan(x[0] || isinf(x[0])))
		return;

	if (W == 0)
		W = dIWidth();
	int NB     = 1 << (PREC - W);
	double step = pow(2.0, -W);

	dIUpdate1_(fold, W, NB, step, x, c, ldx, y);
}

void zIUpdates1(int fold, int W, double complex* x, double complex* c, int ldx, double y) {
	dIUpdate1(fold, W, (double*)x, (double*)(c), 2 * ldx, fabs(y));
	dIUpdate1(fold, W, ((double*)x) + 1, (double*)(c) + 1, 2 * ldx, fabs(y));
}

void zIUpdate1(int fold, int W,double complex* x, double complex* c, int ldx,double complex y) {
	double* tmp = (double*)&y;
	dIUpdate1(fold, W, (double*)x, (double*)(c), 2 * ldx, fabs(tmp[0]));
	dIUpdate1(fold, W, ((double*)x) + 1, (double*)(c) + 1, 2 * ldx, fabs(tmp[1]));
}

