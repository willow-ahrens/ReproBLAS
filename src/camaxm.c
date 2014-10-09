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

float complex camaxm(
	int n,
	float complex* x, int incx,
	float complex* y, int incy
) {
/*
 * Purpose
 *  Compute the maximum absolute value of the element-wise
 *  product of two arrays of single complex nb
 * 
 * Arguments
 *  n     [input] size of input vectors
 *  x     [input] pointer of the first element of the first array
 *  incx  [input] increment between two consecutive elements of X
 *  y     [input] pointer of the first element of the second array
 *  incx  [input] increment between two consecutive elements of Y
 * Return
 *  a single complex floating point number R where
 *   Real(R)      = max_i (abs(Real(x[i]) * Real(y[i])))
 *   Imaginary(R) = max_i (abs(Imaginary(x[i]) * Imaginary(y[i])))
 *
 */
	int i;
	float S0 = 0.0;
	float S1 = 0.0;
	float x0, x1, y0, y1;
	float v0, v1;
	float real, imag;
	float complex ret;

	real = imag = 0.0;

	i = 0;

#ifdef __SSE2__
	float tmp[4] __attribute__ ((aligned(16)));
	__m128 mS, mS2, mv, mv2, my, my2, mAbsMask;
	__m128 mReal, mImag;

	SSE_ABS_MASKS(mAbsMask);

	mReal = _mm_setzero_ps();
	mImag = _mm_setzero_ps();

	for (; i < n; i++, x += incx, y += incy) {
		mv = _mm_loadu_ps((float*)x);
		my = _mm_loadu_ps((float*)y);
		mv2 = _mm_shuffle_ps(mv, mv, _MM_SHUFFLE(0,1,0,1));

		mv  = _mm_mul_ps(mv, my);
		mv2 = _mm_mul_ps(mv2, my);

		mv = _mm_and_ps(mAbsMask, mv);
		mv2 = _mm_and_ps(mAbsMask, mv2);

		mReal = _mm_max_ps(mReal, mv);
		mImag = _mm_max_ps(mImag, mv2);
	}

	_mm_store_ps(tmp, mReal);
	real = tmp[0] > tmp[1] ? tmp[0] : tmp[1];

	_mm_store_ps(tmp, mImag);
	imag = tmp[0] > tmp[1] ? tmp[0] : tmp[1];

#else
	float* ptrx = (float*) x;
	float* ptry = (float*) y;
	incx *= 2;
	incy *= 2;
	for (; i < n; i++, ptrx += incx, ptry += incy) {
		v0 = fabs(ptrx[0] * ptry[0]);
		v1 = fabs(ptrx[1] * ptry[1]);
		real = real > v0 ? real:v0;
		real = real > v1 ? real:v1;

		v0 = fabs(ptrx[1] * ptry[0]);
		v1 = fabs(ptrx[0] * ptry[1]);
		imag = imag > v0 ? imag:v0;
		imag = imag > v1 ? imag:v1;
	}
#endif

	CSET_(ret, real, imag);
	return ret;
}


