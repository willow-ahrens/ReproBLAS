/*
 *  Created   13/10/25   H.D. Nguyen
 */

#ifndef __SSE2__
#include "sIBlas1_SISD.c"
#else

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <xmmintrin.h>

#include "config.h"
#include "Common/Common.h"
#include "sIBlas1.h"

/*
 * SUMMATION IN INDEXEDFP FORMAT
 * REFERENCE IMPLEMENTATION
 */

#define SSUMI2_UNROLL_STEP 4
#define MANUAL_UNROLL

#if defined( SASUMI2 )
void sasumI2(int n, float* v, int incv, int fold, float* sum) {
#elif defined( SNRM2I2 )
void snrm2I2(int n, float* v, int incv, float scale, int fold, float* sum) {
#elif defined( SSUMI2 )
void ssumI2(int n, float* v, int incv, int fold, float* sum) {
#elif defined( SDOTI2 )
void sdotI2(int n, float* v, int incv, float* y, int incy, int fold, float* sum) {
#endif
	
	int i, j;
	__m128 mBUF[MAX_FOLD];
	float tmp[4] __attribute__((aligned(16)));
	__m128 mv0, mv1;
	__m128 mq0, mq1;
	__m128 mS0, mR;
	__m128 mBLP;

	SET_DAZ_FLAG

	// ABSOLUTE MASK
#if defined( SASUMI2 )
	__m128 mAbsMask;
	SSE_ABS_MASKS(mAbsMask);
#elif defined( SNRM2I2 )
	__m128 mScale;
	mScale = _mm_set1_ps(scale);
#elif defined( SDOTI2 )
	__m128 my0, my1;
#endif

	// SET LAST BIT MASK: mBLP = 0.000...1 * 2^0
	SSE_BLP_MASKS(mBLP);

	// EXPAND INITIAL SUM TO BUFFER
	for (j = 0; j < fold; j++) {
		mBUF[j] = _mm_load1_ps(sum + j);
	}

	i = 0;

#ifndef SDOTI2
	if (incv == 1) {
#else
	if (incv == 1 && incy == 1) {
#endif
	while (IS_UNALIGNED(v, 16)) {
		mv0 = _mm_load_ss(v);
#if   defined( SASUMI2 )
		mv0 = _mm_and_ps(mv0, mAbsMask);
#elif defined( SNRM2I2 )
		mv0 = _mm_mul_ss(mv0, mScale);
		mv0 = _mm_mul_ss(mv0, mv0);
#elif defined( SDOTI2 )
		my0 = _mm_load_ss(y);
		mv0 = _mm_mul_ps(mv0, my0);
        y++;
#endif

		for (j = 0; j < fold - 1; j++) {
			mq1     = _mm_or_ps(mv0, mBLP);
			mR      = mBUF[j];
			mBUF[j] = _mm_add_ps(mBUF[j], mq1);
			mR      = _mm_sub_ps(mR, mBUF[j]);
			mv0     = _mm_add_ps(mv0, mR);
		}

		mv0     = _mm_or_ps(mv0, mBLP);
		mBUF[j] = _mm_add_ps(mBUF[j], mv0);
		i ++;
		v++;
	}


#if (SSUMI2_UNROLL_STEP >= 4)
	#ifdef MANUAL_UNROLL
	if (fold == 3) {
		__m128 mT0, mT1, mT2;
		mT0 = mBUF[0];
		mT1 = mBUF[1];
		mT2 = mBUF[2];
		// TODO: check if Y is well aligned
	for (; i < n-7; i+=8, v += 8) {
		mv0 = _mm_load_ps(v);
		mv1 = _mm_load_ps(v + 4);

#if defined( SASUMI2 )
		mv0 = _mm_and_ps(mv0, mAbsMask);
		mv1 = _mm_and_ps(mv1, mAbsMask);
#elif defined( SNRM2I2 )
		mv0 = _mm_mul_ps(mv0, mScale);
		mv0 = _mm_mul_ps(mv0, mv0);
		mv1 = _mm_mul_ps(mv1, mScale);
		mv1 = _mm_mul_ps(mv1, mv1);
#elif defined( SDOTI2 )
		my0 = _mm_load_ps(y);
		my1 = _mm_load_ps(y + 4);
		mv0 = _mm_mul_ps(mv0, my0);
		mv1 = _mm_mul_ps(mv1, my1);
		y += 8;
#endif
		// -- 0 --
			mq0     = _mm_or_ps(mv0, mBLP);
			mq1     = _mm_or_ps(mv1, mBLP);

			mq0     = _mm_add_ps(mT0, mq0);
			mS0     = mT0;
			mq1     = _mm_add_ps(mq0, mq1);

			mS0     = _mm_sub_ps(mS0, mq0);
			mT0     = mq1;
			mq0     = _mm_sub_ps(mq0, mq1);

			mv0     = _mm_add_ps(mv0, mS0);
			mv1     = _mm_add_ps(mv1, mq0);

		// -- 1 --
			mq0     = _mm_or_ps(mv0, mBLP);
			mq1     = _mm_or_ps(mv1, mBLP);

			mq0     = _mm_add_ps(mT1, mq0);
			mS0     = mT1;
			mq1     = _mm_add_ps(mq0, mq1);

			mS0     = _mm_sub_ps(mS0, mq0);
			mT1     = mq1;
			mq0     = _mm_sub_ps(mq0, mq1);

			mv0     = _mm_add_ps(mv0, mS0);
			mv1     = _mm_add_ps(mv1, mq0);

		// -- 2 --
		mv0     = _mm_or_ps(mv0, mBLP);
		mv1     = _mm_or_ps(mv1, mBLP);
		mT2     = _mm_add_ps(mT2, mv0);
		mT2     = _mm_add_ps(mT2, mv1);
	}
		mBUF[0] = mT0;
		mBUF[1] = mT1;
		mBUF[2] = mT2;
	}
	else
	#endif
	for (; i < n-7; i+=8, v += 8) {
		mv0 = _mm_load_ps(v);
		mv1 = _mm_load_ps(v + 4);

#if defined( SASUMI2 )
		mv0 = _mm_and_ps(mv0, mAbsMask);
		mv1 = _mm_and_ps(mv1, mAbsMask);
#elif defined( SNRM2I2 )
		mv0 = _mm_mul_ps(mv0, mScale);
		mv0 = _mm_mul_ps(mv0, mv0);
		mv1 = _mm_mul_ps(mv1, mScale);
		mv1 = _mm_mul_ps(mv1, mv1);
#elif defined( SDOTI2 )
		my0 = _mm_load_ps(y);
		my1 = _mm_load_ps(y + 4);
		mv0 = _mm_mul_ps(mv0, my0);
		mv1 = _mm_mul_ps(mv1, my1);
		y += 8;
#endif

		for (j = 0; j < fold - 1; j++) {
			mq0     = _mm_or_ps(mv0, mBLP);
			mq1     = _mm_or_ps(mv1, mBLP);

			mS0     = mBUF[j];
			mR     = _mm_add_ps(mS0, mq0);
			mBUF[j] = _mm_add_ps(mR, mq1);

			mS0     = _mm_sub_ps(mS0, mR);
			mR     = _mm_sub_ps(mR, mBUF[j]);

			mv0     = _mm_add_ps(mv0, mS0);
			mv1     = _mm_add_ps(mv1, mR);
		}

		mv0     = _mm_or_ps(mv0, mBLP);
		mv1     = _mm_or_ps(mv1, mBLP);
		mBUF[j] = _mm_add_ps(mBUF[j], mv0);
		mBUF[j] = _mm_add_ps(mBUF[j], mv1);
	}
	
#endif

#if (SSUMI2_UNROLL_STEP >= 2)
	for (; i < n-3; i+=4, v += 4) {
		mv0 = _mm_load_ps(v);
#		if defined( SASUMI2 )
		mv0 = _mm_and_ps(mv0, mAbsMask);
#		elif defined( SNRM2I2 )
		mv0 = _mm_mul_ps(mv0, mScale);
		mv0 = _mm_mul_ps(mv0, mv0);
#		elif defined( SDOTI2 )
		my0 = _mm_load_ps(y);
		mv0 = _mm_mul_ps(mv0, my0);
		y += 4;
#endif

		for (j = 0; j < fold - 1; j++) {
			mq1     = _mm_or_ps(mv0, mBLP);
			mR     = mBUF[j];
			mBUF[j] = _mm_add_ps(mBUF[j], mq1);
			mR     = _mm_sub_ps(mR, mBUF[j]);
			mv0     = _mm_add_ps(mv0, mR);
		}

		mv0     = _mm_or_ps(mv0, mBLP);
		mBUF[j] = _mm_add_ps(mBUF[j], mv0);
	}
#endif
	}
	else {
	for (; i < n-3; i+=4, v += 4 * incv) {
		mv0 = _mm_set_ps(v[0], v[incv], v[2*incv], v[3*incv]);

#		if defined( SASUMI2 )
		mv0 = _mm_and_ps(mv0, mAbsMask);
#		elif defined( SNRM2I2 )
		mv0 = _mm_mul_ps(mv0, mScale);
		mv0 = _mm_mul_ps(mv0, mv0);
#		elif defined( SDOTI2 )
		my0 = _mm_set_ps(y[0], y[incy], y[2*incy], y[3*incy]);
		mv0 = _mm_mul_ps(mv0, my0);
		y += 4 * incy;
#		endif

		for (j = 0; j < fold - 1; j++) {
			mq1     = _mm_or_ps(mv0, mBLP);
			mR     = mBUF[j];
			mBUF[j] = _mm_add_ps(mBUF[j], mq1);
			mR     = _mm_sub_ps(mR, mBUF[j]);
			mv0     = _mm_add_ps(mv0, mR);
		}

		mv0     = _mm_or_ps(mv0, mBLP);
		mBUF[j] = _mm_add_ps(mBUF[j], mv0);
	}
	}
	// last one
	for (; i < n; i++, v+=incv) {
		mv0 = _mm_load_ss(v);
#		if defined( SASUMI2 )
		mv0 = _mm_and_ps(mv0, mAbsMask);
#		elif defined( SNRM2I2 )
		mv0 = _mm_mul_ss(mv0, mScale);
		mv0 = _mm_mul_ss(mv0, mv0);
#		elif defined( SDOTI2 )
		my0 = _mm_load_ss(y);
		mv0 = _mm_mul_ss(mv0, my0);
		y += incy;
#		endif

		for (j = 0; j < fold - 1; j++) {
			mq1     = _mm_or_ps(mv0, mBLP);
			mq1 = mv0;
			mR     = mBUF[j];
			mBUF[j] = _mm_add_ps(mBUF[j], mq1);
			mR     = _mm_sub_ps(mR, mBUF[j]);
			mv0     = _mm_add_ps(mv0, mR);
		}

		mv0     = _mm_or_ps(mv0, mBLP);
		mBUF[j] = _mm_add_ps(mBUF[j], mv0);
	}


	// TODO SSE3 horizontal add _mm_hadd_ps
	for (j = 0; j < fold; j++) {
		mR  = _mm_load1_ps(sum + j);
		mq1 = _mm_load_ss(sum + j);
		mR  = _mm_sub_ps(mR, mq1);
		mBUF[j] = _mm_sub_ps(mBUF[j], mR);
#ifdef ___SSE3__
		mBUF[j] = _mm_hadd_ps(mBUF[j], mBUF[j]);
		mBUF[j] = _mm_hadd_ps(mBUF[j], mBUF[j]);
		_mm_store_ss(sum + j, mBUF[j]);
#else
		_mm_store_ps(tmp, mBUF[j]);
		sum[j] = tmp[0] + tmp[1] + tmp[2] + tmp[3];
#endif
	}
	RESET_DAZ_FLAG
}
#endif // ifndef __SSE2__
