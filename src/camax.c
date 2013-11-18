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

#define CAMAX_UNROLL  2

float complex camax(int n, float complex* vd, int inc) {
/*
 * Purpose
 *  Compute the maximum absolute value of an array of single complex nb
 * 
 * Arguments
 *  n     [input] size of input vector
 *  vd    [input] pointer of the first element of the array
 *  inc   [input] increment between two consecutive elements
 * Return
 *  a single complex floating point number R where
 *   Real(R)      = max_i (abs(Real(v[i])))
 *   Imaginary(R) = max_i (abs(Imaginary(v[i])))
 *
 */
	int i, j;
	float S0 = 0.0;
	float S1 = 0.0;
	float v0, v1;
	float tmp[4] __attribute__ ((aligned(16)));
	float complex ret;
	float* ptr = (float*) &ret;

	i = 0;

#ifdef __SSE2__
	__m128 mS, mv, mAbsMask, mS1, mv1;
	mS  = _mm_setzero_ps();
	mS1 = _mm_setzero_ps();
	SIMD_ABS_MASKS(mAbsMask);
#endif
	
	// convert to float pointer
	float* v = (float*) vd;
	int pad = 0;
	if (inc == 1) {
#ifdef __SSE2__
		if (IS_UNALIGNED(v)) {
			S0 = fabs(v[0]);
			v += inc;
			i ++;
			if (IS_UNALIGNED(v)) {
				S1 = fabs(v[0]);
				v += inc;
				i ++;
				if (IS_UNALIGNED(v)) {
					v0 = fabs(v[0]);
					S0 = S0 < v0 ? v0 : S0;
					v += inc;
					i ++;
				}
			}
		}
		pad = i;


#if (CAMAX_UNROLL >= 4)
		for (; i < 2*n-7; i+=8) {
			mv  = _mm_load_ps(v);
			mv1 = _mm_load_ps(v+4);

			mv  = _mm_and_ps(mAbsMask, mv);
			mv1 = _mm_and_ps(mAbsMask, mv1);

			mS  = _mm_max_ps(mS, mv);
			mS1 = _mm_max_ps(mS1, mv1);
			v += 8;
		}

		mS = _mm_max_ps(mS, mS1);
#endif
		for (; i < 2*n-3; i+=4,v+=4) {
			mv = _mm_load_ps(v);
			mv = _mm_and_ps(mAbsMask, mv);

			mS = _mm_max_ps(mS, mv);
		}

		_mm_store_ps(tmp, mS);

		for (j = 0; i < 2*n; i++, v++, j++) {
			v0 = fabs(v[0]);
			tmp[j] = tmp[j] < v0 ? v0 : tmp[j];
		}

		tmp[0] = tmp[0] < tmp[2] ? tmp[2] : tmp[0];
		tmp[1] = tmp[1] < tmp[3] ? tmp[3] : tmp[1];

		if (pad) {
			S0 = S0 < tmp[1] ? tmp[1] : S0;
			S1 = S1 < tmp[0] ? tmp[0] : S1;
		}
		else {
			S0 = S0 < tmp[0] ? tmp[0] : S0;
			S1 = S1 < tmp[1] ? tmp[1] : S1;
		}
		ptr[0] = S0;
		ptr[1] = S1;
		return ret;
#else
		float S2 = 0.0, S3 = 0.0;
		float v2, v3;
		for (; i < 2*n-3; i+=4, v+=4) {
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

		if (i < 2*n-1) {
			v0 = fabs(v[0]);
			v1 = fabs(v[1]);

			S0 = MAX(S0, v0);
			S1 = MAX(S1, v1);
			i+=2; v+=2;
		}

		if (i < 2*n) {
			v0 = fabs(v[0]);
			S0 = MAX(S0, v0);

			i++; v++;
		}
		// IF DATA IS UNALIGNED, NEED TO SHUFFLE
		if (pad % 2 == 1) {
			v0 = fabs(v[0]);
			S0 = MAX(S0, v0);

			ptr[1] = S0;
			ptr[0] = S1;
		}
		else {
			ptr[0] = S0;
			ptr[1] = S1;
		}
		return ret;
#endif
	}

	// DATA ARE NOT CONTIGUOUS

	inc *= 2;
#ifdef __SSE2__
	for (; i < n; i++,v+=inc) {
		mv = _mm_loadu_ps(v);
		mv = _mm_and_ps(mAbsMask, mv);

		mS = _mm_max_ps(mS, mv);
	}

	_mm_store_ps(tmp, mS);
	S0 = tmp[0];
	S1 = tmp[1];
#else
	for (; i < n; i++, v += inc) {
		v0 = fabs(v[0]);
		v1 = fabs(v[1]);

		S0 = max(S0, v0);
		S1 = max(S1, v1);
	}
#endif
	
	ptr[0] = S0;
	ptr[1] = S1;

	return ret;
}

