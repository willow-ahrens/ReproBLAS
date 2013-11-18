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

#define DAMAX_UNROLL  2

double complex zamax(int n, double complex* vd, int inc) {
	int i;
	double S0 = 0.0;
	double S1 = 0.0;
	double v0, v1;
	double complex ret[1] __attribute__ ((aligned(16)));
	double* tmp = (double*)ret;

	i = 0;

#ifdef __SSE2__
	__m128d mS, mv, mAbsMask, mS1, mv1;
	mS  = _mm_setzero_pd();
	mS1 = _mm_setzero_pd();
	SIMD_ABS_MASKD(mAbsMask);
#endif
	
	int unaligned = IS_UNALIGNED(vd);
	double* v = (double*) vd;
	if (inc == 1) {
#ifdef __SSE2__
		// convert to double pointer
		if (unaligned) {
			mS = _mm_loadh_pd(_mm_setzero_pd(), v);
//			S1 = fabs(v[0]);
			v ++;
			i = 1;
		}


#if (DAMAX_UNROLL >= 4)
		for (; i < 2*n-3; i+=4, v+=4) {
			mv  = _mm_load_pd(v);
			mv1 = _mm_load_pd(v+2);

			mv  = _mm_and_pd(mAbsMask, mv);
			mv1 = _mm_and_pd(mAbsMask, mv1);

			mS  = _mm_max_pd(mS, mv);
			mS1 = _mm_max_pd(mS1, mv1);
		}

		mS = _mm_max_pd(mS, mS1);
#endif
		for (; i < 2*n-1; i+=2,v+=2) {
			mv = _mm_load_pd(v);
			mv = _mm_and_pd(mAbsMask, mv);

			mS = _mm_max_pd(mS, mv);
		}

		if (i < 2*n) {
			mv = _mm_load_sd(v);
			mv = _mm_and_pd(mAbsMask, mv);

			mS = _mm_max_pd(mS, mv);
		}

		if (unaligned) {
			// SHUFFLE
			mS = _mm_shuffle_pd(mS, mS, _MM_SHUFFLE2(0,1));
		}

		_mm_store_pd(tmp, mS);
		return ret[0];
#else
#if DAMAX_UNROLL >= 4
		double S2 = 0.0, S3 = 0.0;
		double v2, v3;
		for (; i < n-1; i+=2, v+=4) {
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
		for (; i < n; i++, v+=2) {
			v0 = fabs(v[0]);
			v1 = fabs(v[1]);

			S0 = max(S0, v0);
			S1 = max(S1, v1);
		}
		tmp[0] = S0;
		tmp[1] = S1;
		return ret[0];
#endif
	}

	// DATA ARE NOT CONTIGUOUS

	inc = inc * 2;
#ifdef __SSE2__
	if (unaligned) {
#if (DAMAX_UNROLL >= 4)
		for (; i < n-1; i+=2, v += 2*inc) {
			mv  = _mm_loadu_pd(v);
			mv1 = _mm_loadu_pd(v+inc);

			mv  = _mm_and_pd(mAbsMask, mv);
			mv1 = _mm_and_pd(mAbsMask, mv1);

			mS  = _mm_max_pd(mS, mv);
			mS1 = _mm_max_pd(mS1, mv1);
		}

		mS = _mm_max_pd(mS, mS1);
#endif
		for (; i < n; i++,v+=inc) {
			mv = _mm_loadu_pd(v);
			mv = _mm_and_pd(mAbsMask, mv);

			mS = _mm_max_pd(mS, mv);
		}
	}
	else {
#if (DAMAX_UNROLL >= 4)
		for (; i < n-1; i+=2, v += 2*inc) {
			mv  = _mm_load_pd(v);
			mv1 = _mm_load_pd(v+inc);

			mv  = _mm_and_pd(mAbsMask, mv);
			mv1 = _mm_and_pd(mAbsMask, mv1);

			mS  = _mm_max_pd(mS, mv);
			mS1 = _mm_max_pd(mS1, mv1);
		}

		mS = _mm_max_pd(mS, mS1);
#endif
		for (; i < n; i++,v+=inc) {
			mv = _mm_load_pd(v);
			mv = _mm_and_pd(mAbsMask, mv);

			mS = _mm_max_pd(mS, mv);
		}
	}
	_mm_store_pd(tmp, mS);
	return ret[0];
#else
	for (; i < n; i++, v += inc) {
		v0 = fabs(v[0]);
		v1 = fabs(v[1]);

		S0 = max(S0, v0);
		S1 = max(S1, v1);
	}
	tmp[0] = S0;
	tmp[1] = S1;
	return ret[0];
#endif
}

