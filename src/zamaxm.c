/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <float.h>
#include <math.h>
#ifdef __SSE2__
#	include <emmintrin.h>
#endif
#include "Common/Common.h"
#include "rblas1.h"
#include "types.h"

double complex zamaxm(int n, double complex* x, int incx, double complex* y, int incy) {
	int i;
	double S0 = 0.0;
	double S1 = 0.0;
	double x0, x1, y0, y1;
	double v0, v1;
	double real, imag;
	double complex ret;

	real = imag = 0.0;

	i = 0;

#ifdef __SSE2__
	double tmp[2] __attribute__ ((aligned(16)));
	__m128d mS, mS2, mv, mv2, my, my2, mAbsMask;
	__m128d mReal, mImag;

	mv = _mm_set1_pd(1);

	SIMD_ABS_MASKD(mAbsMask);

	mReal = _mm_setzero_pd();
	mImag = _mm_setzero_pd();

	for (; i < n; i++, x += incx, y += incy) {
		mv = _mm_loadu_pd((double*)x);
		my = _mm_loadu_pd((double*)y);
		mv2 = _mm_loadr_pd((double*)x);

		mv  = _mm_mul_pd(mv, my);
		mv2 = _mm_mul_pd(mv2, my);

		mv = _mm_and_pd(mAbsMask, mv);
		mv2 = _mm_and_pd(mAbsMask, mv2);

		mReal = _mm_max_pd(mReal, mv);
		mImag = _mm_max_pd(mImag, mv2);
	}

	_mm_store_pd(tmp, mReal);
	real = tmp[0] > tmp[1] ? tmp[0] : tmp[1];

	_mm_store_pd(tmp, mImag);
	imag = tmp[0] > tmp[1] ? tmp[0] : tmp[1];

#else
	for (; i < n; i++, x += incx, y+= incy) {
		v0 = fabs(x[0].real * y[0].real);
		v1 = fabs(x[0].imag * y[0].imag);
		real = real > v0 ? real:v0;
		real = real > v1 ? real:v1;

		v0 = fabs(x[0].real * y[0].imag);
		v1 = fabs(x[0].imag * y[0].real);
		imag = imag > v0 ? imag:v0;
		imag = imag > v1 ? imag:v1;
	}
#endif

	ZSET_(ret, real, imag);
	return ret;
}


