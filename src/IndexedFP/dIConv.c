/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

void dconv2I1(int fold, int W, double x, double* rep, double* C, int inc) {
	int i;
	double q;
	if (x == 0.0) {
		for (i = 0; i < fold; i++) {
			rep[i*inc] = 0.0;
			C  [i*inc] = 0.0;
		}
		return;
	}
	dIBoundary(fold, W, fabs(x), rep, inc);
	double M;
	for (i = 0; i < fold; i++, rep += inc, C += inc) {
		// high-order part
		C[0] = 0.0;

		M = rep[0];
		q = M + x;
		rep[0] = q;
		q -= M;
		x -= q;
	}
}

I_double dconv2I(double x) {
	I_double ret;
	dconv2I1(DEFAULT_FOLD, 0, x, ret.m, ret.c, 1);
	return ret;
}

double Iconv2d1(int fold, double* rep, double* carry, int inc) {
	int i;
	double ret = 0.0;

	// CHECK FOR NAN OR INFINITY
	if (isinf(rep[0]) || isnan(rep[0]))
		return rep[0];

	if (rep[0] == 0.0) {
		return 0.0;
	}


	double M;
	
	// TODO: SCALING TO AVOID OVERFLOW
	fold = (fold < 1) ? 1 : fold;

#if 1
	for (i = 0; i < fold; i++, rep += inc, carry += inc) {
		M = ufp(rep[0]);
		ret += (rep[0] + (carry[0] - 6) * M * 0.25);
	}
#else // USING DOUBLE-DOUBLE
	double rdd[4];
	rdd[0] = rdd[1] = rdd[2] = rdd[3] = 0.0;
	for (i = 0; i < fold; i++) {
		M = ufp(rep[i]) * 0.25;
		ddpd(rdd, (rep[i + fold] - 6) * M);
		ddpd(rdd, rep[i]);
	}
	ret = rdd[0] + rdd[1] + rdd[2] + rdd[3];
#endif

	return ret;
}

double complex Iconv2z1(int fold, double complex* rep, double complex* C, int inc) {
	double complex ret;
	double* rptr = (double*) &ret;
	inc *= 2;
	rptr[0] = Iconv2d1(fold, (double*)rep, (double*)C, inc);
	rptr[1] = Iconv2d1(fold, (double*)rep + 1, (double*)C + 1, inc);
	return ret;
}

void zconv2I1(int fold, int W, double complex X, double complex* rep, double complex* C, int inc) {
	double* xx = (double*) &X;
	double* dptr = (double*) rep;
	double* cptr = (double*) C;
	inc *= 2;
	dconv2I1(fold, W, xx[0], dptr, cptr, inc);
	dconv2I1(fold, W, xx[1], dptr+1, cptr + 1, inc);
}

I_double_Complex zconv2I(double complex x) {
	I_double_Complex ret;
	zconv2I1(DEFAULT_FOLD, 0, x, (double complex*)ret.m, (double complex*) ret.c, 1);
	return ret;
}

