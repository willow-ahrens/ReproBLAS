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

void fconv2I1(int fold, float x, float* mantissa, float* carry, int inc) {
	int i;
	float q;
	for (i = 0; i < fold; i++) carry[i*inc] = 0;
	if (x == 0.0) {
		for (i = 0; i < fold; i++)
			mantissa[i*inc] = 0.0;
		return;
	}
	sIBoundary_(fold, fabs(x), mantissa, inc);
	float M;
	for (i = 0; i < fold; i++, mantissa += inc) {
		M = mantissa[0];
		q = M + x;
		mantissa[0] = q;
		q -= M;
		x -= q;
	}
}

I_float fconv2I(float x) {
	I_float ret;
	fconv2I1(DEFAULT_FOLD, x, ret.m, ret.c, 1);
	return ret;
}

float Iconv2f1(int fold, float* mantissa, float* carry, int inc) {
	int i;

	// CHECK FOR NAN OR INFINITY
	if (isinf(mantissa[0]) || isnan(mantissa[0]))
		return mantissa[0];

	if (mantissa[0] == 0.0) {
		return 0.0;
	}

	float M;
	float step;
	double ret = 0.0;

#ifdef DEBUG
	fprintf(stdout, "CONVERT BACK [{%d}", fold);
	for (i = 0; i < fold; i++)
		fprintf(stdout, "{%g,%g}", carry[i], mantissa[i]);
	fprintf(stdout, "] ");
#endif

	// TODO: SCALING TO AVOID OVERFLOW

	fold = (fold < 1) ? 1 : fold;

#if 1
	for (i = 0; i < fold; i++, mantissa += inc, carry += inc) {
		M = ufpf(mantissa[0]);
		ret += (double)(carry[0] - 6) * (double)(M * 0.25);
		ret += mantissa[0];
	}
#else // USING DOUBLE-DOUBLE
		float rdd[4];
		rdd[0] = rdd[1] = rdd[2] = rdd[3] = 0.0;
		for (i = 1; i < fold + 1; i++) {
			ddpd(rdd, mantissa[i + fold] * M);
			ddpd(rdd, mantissa[i]);
			M *= step;
		}
		ret = rdd[0] + rdd[1] + rdd[2] + rdd[3];
#endif

#ifdef DEBUG
	fprintf(stdout, " -> %g \n", ret);
#endif
	return (float)ret;
}

float complex Iconv2c1(int fold, float complex* mantissa, float* C, int inc) {
	float complex ret;
	float* rptr = (float*) &ret;
	inc *= 2;
	rptr[0] = Iconv2f1(fold, (float*)mantissa, (float*)C, inc);
	rptr[1] = Iconv2f1(fold, (float*)mantissa + 1, (float*)C + 1, inc);
	return ret;
}

void cconv2I1(int fold, float complex x, float complex* m, float* c, int inc) {
	float* xptr = (float*) &(x);
	fconv2I1(fold, xptr[0], (float*)m  , c  , 2*inc);
	fconv2I1(fold, xptr[1], (float*)m+1, c+1, 2*inc);
}

I_float_Complex cconv2I(float complex x) {
	I_float_Complex ret;
	cconv2I1(DEFAULT_FOLD, x, (float complex*)ret.m, ret.c, 1);
	return ret;
}
