/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <float.h>
#include <math.h>
#ifdef __SSE__
#	include <xmmintrin.h>
#endif
#include "Common/Common.h"

#define SAMAX_UNROLL 8 

float samax(int n, float* v, int inc) {
	int i;
	float S0 = 0.0;
	float S1 = 0.0;
	float v0, v1;
	float tmp[4] __attribute__ ((aligned(16)));


	i = 0;

	
	if (inc == 1) {
		while (IS_UNALIGNED(v)) {
			v0 = fabs(v[0]);
			if (S1 < v0) S1 = v0;
			v += inc;
			i ++;
		}
#ifdef __SSE__
		__m128 mS, mv, mAbsMask, mS1, mv1;

		SIMD_ABS_MASKS(mAbsMask);

		mS  = _mm_setzero_ps();
		mS1 = _mm_setzero_ps();

#if (SAMAX_UNROLL >= 8)
		for (; i < n-7; i+=8, v += 8) {
			mv  = _mm_load_ps(v);
			mv1 = _mm_load_ps(v+4);

			mv  = _mm_and_ps(mAbsMask, mv);
			mv1 = _mm_and_ps(mAbsMask, mv1);

			mS  = _mm_max_ps(mS, mv);
			mS1 = _mm_max_ps(mS1, mv1);
		}

		mS = _mm_max_ps(mS, mS1);
#else
		for (; i < n-3; i+=4,v+=4) {
			mv = _mm_load_ps(v);
			mv = _mm_and_ps(mAbsMask, mv);

			mS = _mm_max_ps(mS, mv);
		}
#endif

		_mm_store_ps(tmp, mS);
		S0 = MAX(tmp[0], tmp[1]);
		S0 = MAX(S0, S1);
		S1 = MAX(tmp[2], tmp[3]);
		S0 = MAX(S0, S1);
#else
#if SAMAX_UNROLL >= 4
		float S2 = 0.0, S3 = 0.0;
		float v2, v3;
		for (; i < n-3; i+=4, v+=4) {
			v0 = fabs(v[0]);
			v1 = fabs(v[1]);
			v2 = fabs(v[2]);
			v3 = fabs(v[3]);

			S0 = MAX(S0, v0);
			S1 = MAX(S1, v1);
			S2 = MAX(S2, v2);
			S3 = MAX(S3, v3);
		}
		S0 = MAX(S0, S2);
		S1 = MAX(S1, S3);
#endif
		/*
		for (; i < n-1; i+=2, v+=2) {
			v0 = fabs(v[0]);
			v1 = fabs(v[1]);

			S0 = MAX(S0, v0);
			S1 = MAX(S1, v1);
		}
		S0 = MAX(S0, S1);
		*/
#endif
	}

	for (; i < n; i++, v += inc) {
		v0 = fabs(v[0]);
		S0 = MAX(S0, v0);
	}
	return S0;
}

