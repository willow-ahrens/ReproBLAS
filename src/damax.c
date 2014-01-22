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

#define DAMAX_UNROLL 4 

double damax(int n, double* v, int inc) {
	int i;
	double S0 = 0.0;
	double S1 = 0.0;
	double v0, v1;
	double tmp[2] __attribute__ ((aligned(16)));


	i = 0;

	
	if (inc == 1) {
		if (IS_UNALIGNED(v, 16)) {
			S1 = fabs(v[0]);
			v += inc;
			i = 1;
		}
#ifdef __SSE2__
		__m128d mS, mv, mAbsMask, mS1, mv1;

		mv = _mm_set1_pd(1);
		mAbsMask = _mm_set1_pd(-1);
		mAbsMask = _mm_xor_pd(mAbsMask, mv);
		mS = _mm_cmpeq_pd(mv, mv);
		mAbsMask = _mm_xor_pd(mAbsMask, mS);
		mS  = _mm_setzero_pd();
		mS1 = _mm_setzero_pd();

#if (DAMAX_UNROLL >= 4)
		for (; i < n-3; i+=4) {
			mv  = _mm_load_pd(v);
			mv1 = _mm_load_pd(v+2);

			mv  = _mm_and_pd(mAbsMask, mv);
			mv1 = _mm_and_pd(mAbsMask, mv1);

			mS  = _mm_max_pd(mS, mv);
			mS1 = _mm_max_pd(mS1, mv1);
			v += 4;
		}

		mS = _mm_max_pd(mS, mS1);
#else
		for (; i < n-1; i+=2,v+=2) {
			mv = _mm_load_pd(v);
			mv = _mm_and_pd(mAbsMask, mv);

			mS = _mm_max_pd(mS, mv);
		}
#endif

		_mm_store_pd(tmp, mS);
		S0 = MAX(tmp[0], tmp[1]);
		S0 = MAX(S0, S1);
#else
#if DAMAX_UNROLL >= 4
		double S2 = 0.0, S3 = 0.0;
		double v2, v3;
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
		for (; i < n-1; i+=2, v+=2) {
			v0 = fabs(v[0]);
			v1 = fabs(v[1]);

			S0 = MAX(S0, v0);
			S1 = MAX(S1, v1);
		}
		S0 = MAX(S0, S1);
#endif
	}

	for (; i < n; i++, v += inc) {
		v0 = fabs(v[0]);
		S0 = MAX(S0, v0);
	}
	return S0;
}

