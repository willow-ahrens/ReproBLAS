/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <xmmintrin.h>
#include "Common/Common.h"

#define SAMAX_UNROLL 8 

float samaxm(int n, float* x, int incx, float* y, int incy) {
	int i;
	float S0 = 0.0;
	float S1 = 0.0;
	float v0;
	float tmp[4] __attribute__ ((aligned(16)));
	__m128 mS, mS2, mv, mv2, my, my2, mAbsMask;

	SSE_ABS_MASKS(mAbsMask);

	i = 0;

	while (IS_UNALIGNED(y, 16)) {
		v0 = fabs(x[0] * y[0]);
		S0 = S0 > v0 ? S0:v0;
		x++;
		y++;
		i ++;
		x ++;
		y ++;
	}

	if (incx == 1 && incy == 1) {
#ifdef __SSE__	
		mS = _mm_setzero_ps();
		mS2 = _mm_setzero_ps();

#if (SAMAX_UNROLL >= 8)
		for (; i < n-7; i+=8, x += 8, y+=8) {
			mv = _mm_load_ps(x);
			my = _mm_load_ps(y);
			mv2 = _mm_load_ps(x + 4);
			my2 = _mm_load_ps(y + 4);

			mv  = _mm_mul_ps(mv, my);
			mv2 = _mm_mul_ps(mv2, my2);

			mv = _mm_and_ps(mAbsMask, mv);
			mv2 = _mm_and_ps(mAbsMask, mv2);

			mS = _mm_max_ps(mS, mv);
			mS2 = _mm_max_ps(mS2, mv2);
		}
		mS = _mm_max_ps(mS, mS2);
#endif

#if (SAMAX_UNROLL >= 4)
		for (; i < n-3; i+=4, x += 4, y+=4) {
			mv = _mm_load_ps(x);
			my = _mm_load_ps(y);
			mv  = _mm_mul_ps(mv, my);
			mv = _mm_and_ps(mAbsMask, mv);
			mS = _mm_max_ps(mS, mv);
		}
#endif

		_mm_store_ps(tmp, mS);
		S1 = tmp[0] > tmp[1] ? tmp[0] : tmp[1];
		S0 = S0 > S1 ? S0 : S1;
		S1 = tmp[2] > tmp[3] ? tmp[2] : tmp[3];
		S0 = S0 > S1 ? S0 : S1;
#else
		double S2 = 0.0, S3 = 0.0;
		double v1, v2, v3;
		for (; i < n - 3; i+=4, x+=4, y+=4) {
			v0 = fabs(x[0] * y[0]);
			v1 = fabs(x[1] * y[1]);
			v2 = fabs(x[2] * y[2]);
			v3 = fabs(x[3] * y[3]);

			S0 = S0 > v0 ? S0:v0;
			S1 = S1 > v1 ? S1:v1;
			S2 = S2 > v2 ? S2:v2;
			S3 = S3 > v3 ? S3:v3;
		}

		S0 = S0 > S1 ? S0:S1;
		S2 = S2 > S3 ? S2:S3;
		S0 = S0 > S2 ? S0:S2;
#endif
	}
	for (; i < n; i++, x+=incx, y+=incy) {
		v0 = fabs(x[0] * y[0]);
		S0 = S0 > v0 ? S0:v0;
	}
	return S0;
}


