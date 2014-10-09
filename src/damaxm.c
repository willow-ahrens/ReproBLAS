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

double damaxm(int n, double* x, int incx, double* y, int incy) {
	int i;
	double S0 = 0.0;
	double S1 = 0.0;
	double x0, x1, y0, y1;
	double v0, v1;

	i = 0;

	if (incx == 1 && incy == 1) {
#ifdef __SSE2__
	double tmp[2] __attribute__ ((aligned(16)));
	__m128d mS, mS2, mv, mv2, my, my2, mAbsMask;

	mv = _mm_set1_pd(1);

	SSE_ABS_MASKD(mAbsMask);

		if (IS_UNALIGNED(x, 16)) {
			S0 = fabs(x[0] * y[0]);
			i = 1;
			x++; y++;
		}

		mS = _mm_setzero_pd();
		mS2 = _mm_setzero_pd();

		if (IS_ALIGNED(y, 16)) {
		for (; i < n-3; i+=4, x += 4, y += 4) {
			mv = _mm_load_pd(x);
			my = _mm_load_pd(y);
			mv2 = _mm_load_pd(x + 2);
			my2 = _mm_load_pd(y + 2);

			mv  = _mm_mul_pd(mv, my);
			mv2 = _mm_mul_pd(mv2, my2);

			mv = _mm_and_pd(mAbsMask, mv);
			mv2 = _mm_and_pd(mAbsMask, mv2);

			mS = _mm_max_pd(mS, mv);
			mS2 = _mm_max_pd(mS2, mv2);
		}
		}
		else {
		for (; i < n-3; i+=4, x += 4, y += 4) {
			mv = _mm_load_pd(x);
			my = _mm_loadu_pd(y);
			mv2 = _mm_load_pd(x + 2);
			my2 = _mm_loadu_pd(y + 2);

			mv  = _mm_mul_pd(mv, my);
			mv2 = _mm_mul_pd(mv2, my2);

			mv = _mm_and_pd(mAbsMask, mv);
			mv2 = _mm_and_pd(mAbsMask, mv2);

			mS = _mm_max_pd(mS, mv);
			mS2 = _mm_max_pd(mS2, mv2);
		}
		}

		mS = _mm_max_pd(mS, mS2);
		_mm_store_pd(tmp, mS);
		S1 = tmp[0] > tmp[1] ? tmp[0] : tmp[1];
		S0 = S0 > S1 ? S0 : S1;

#else
	for (; i < n-1; i+=2) {
		v0 = fabs(x[0] * y[0]);
		v1 = fabs(x[1] * y[1]);

		S0 = S0 > v0 ? S0:v0;
		S1 = S1 > v1 ? S1:v1;

		x += 2;
		y += 2;
	}
	S0 = S0 > S1 ? S0:S1;
#endif
	}

	for (; i < n; i++) {
		v0 = fabs(x[0] * y[0]);
		S0 = S0 > v0 ? S0:v0;
		x += incx; y += incy;
	}
	return S0;
}


