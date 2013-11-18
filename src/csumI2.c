/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <emmintrin.h>

#include "config.h"
#include "Common/Common.h"
#include "cIBlas1.h"
#include "rblas1.h"

#define MANUAL_UNROLL

#define CSUMI2_UNROLL_STEP 2

/*
 * Reference SIMD implementation for 3-fold binning sum/asum/nrm2/dot
 */

#ifdef MANUAL_UNROLL

#if defined( SCASUMI2 )
void scasumI2_k3(int n, float complex* v, int incv, float complex* sum) {
#endif
#if defined( CSUMI2 )
void csumI2_k3(int n, float complex* v, int incv, float complex* sum) {
#endif
#if defined( SCNRM2I2 )
void scnrm2I2_k3(int n, float complex* v, int incv, float scale, float complex* sum) {
#endif

	int i;

	__m128 mS1, mS2;
	__m128 mT1, mT2;
	__m128 mS1p, mS2p;
	__m128 mL1, mL2;
	__m128 mL1p, mL2p;
	__m128 mv0, mv1;
	float complex tmp[2] __attribute__ ((aligned(16)));

	// GENERATING MASK
	// mBLP = 0.000...1 * 2^0
	__m128 mBLP; // BIT IN THE LAST PLACE MASK
	SIMD_BLP_MASKS(mBLP);

#if defined( SCASUMI2 )
	__m128 mAbsMask;
	SIMD_ABS_MASKS(mAbsMask);
#endif
#ifdef SCNRM2I2
	__m128 mScale;
	mScale = _mm_set1_ps(scale);
#endif

	mS1 = _mm_loadu_ps((float*)sum);
	mS2 = mS1;

	mS1p = mS1;
	mS2p = mS2;

	mT1 = _mm_loadu_ps((float*)(sum + 2));
	mT2 = mT1;

	mL1 = _mm_loadu_ps((float*)(sum + 1));
	mL2 = mL1;

	mL1p = mL1;
	mL2p = mL2;

	i = 0;
	
#if (CSUMI2_UNROLL_STEP >= 4)
	for (; i < n-1; i+=2, v+=2*incv) {
		mv0 = _mm_loadu_ps((float*)v);
		mv1 = _mm_loadu_ps((float*)(v + incv));

#		ifdef SCASUMI2
		// ABSOLUTE VALUE
		mv0 = _mm_and_ps(mv0, mAbsMask);
		mv1 = _mm_and_ps(mv1, mAbsMask);
#		endif
#		ifdef SCNRM2I2
		mv0 = _mm_mul_ps(mv0, mScale);
		mv0 = _mm_mul_ps(mv0, mv0);
		mv1 = _mm_mul_ps(mv1, mScale);
		mv1 = _mm_mul_ps(mv1, mv1);
#endif
		// CHANGE LAST BIT TO 1
		mS1p = _mm_or_ps(mv0, mBLP);
		mS2p = _mm_or_ps(mv1, mBLP);

		//-----
		mS1p  = _mm_add_ps(mS1p, mS1);
		mS1  = _mm_sub_ps(mS1, mS1p);
		mv0  = _mm_add_ps(mv0, mS1);
		mS1  = mS1p;

		mL1p = _mm_or_ps(mv0, mBLP);
		mL1p  = _mm_add_ps(mL1p, mL1);
		mL1  = _mm_sub_ps(mL1, mL1p);
		mv0  = _mm_add_ps(mv0, mL1);
		mv0 = _mm_or_ps(mv0, mBLP);
		mL1  = mL1p;

		mT1  = _mm_add_ps(mT1, mv0);
		//----

		mS2p  = _mm_add_ps(mS2p, mS2);
		mS2  = _mm_sub_ps(mS2, mS2p);
		mv1  = _mm_add_ps(mv1, mS2);
		mS2  = mS2p;

		mL2p = _mm_or_ps(mv1, mBLP);
		mL2p  = _mm_add_ps(mL2p, mL2);
		mL2  = _mm_sub_ps(mL2, mL2p);
		mv1  = _mm_add_ps(mv1, mL2);
		mv1 = _mm_or_ps(mv1, mBLP);
		mL2  = mL2p;

		mT2  = _mm_add_ps(mT2, mv1);
		//-----
	}

	mS1p = _mm_loadu_ps((float*)sum);
	mL1p = _mm_loadu_ps((float*)(sum+1));
	mS2p = _mm_loadu_ps((float*)(sum+2));

	mS2 = _mm_sub_ps(mS2, mS1p);
	mL2 = _mm_sub_ps(mL2, mL1p);
	mT2 = _mm_sub_ps(mT2, mS2p);

	mS1 = _mm_add_ps(mS1, mS2);
	mL1 = _mm_add_ps(mL1, mL2);
	mT1 = _mm_add_ps(mT1, mT2);

	mS1p = mS1;
	mL1p = mL1;
#endif

#if (CSUMI2_UNROLL_STEP >= 2)
#if (CSUMI2_UNROLL_STEP >= 4)
	if (i < n) {
		i += 1;
#else
	for (; i < n; i+=1) {
#endif
		mv0 = _mm_loadu_ps((float*)v);

#if defined( SCASUMI2 )
		mv0 = _mm_and_ps(mv0, mAbsMask);
#endif
#		ifdef SCNRM2I2
		mv0 = _mm_mul_ps(mv0, mScale);
		mv0 = _mm_mul_ps(mv0, mv0);
#endif
		mS1p = _mm_or_ps(mv0, mBLP);

		//-----
		mS1p  = _mm_add_ps(mS1p, mS1);
		mS1  = _mm_sub_ps(mS1, mS1p);
		mv0  = _mm_add_ps(mv0, mS1);
		mS1  = mS1p;

		mL1p = _mm_or_ps(mv0, mBLP);

		mL1p  = _mm_add_ps(mL1p, mL1);
		mL1  = _mm_sub_ps(mL1, mL1p);
		mv0  = _mm_add_ps(mv0, mL1);
		mv0 = _mm_or_ps(mv0, mBLP);
		mL1  = mL1p;

		mT1  = _mm_add_ps(mT1, mv0);
		//----
		v += incv;
	}
#endif

	_mm_store_ps((float*)tmp, mS1);
	sum[0] = tmp[0];

	_mm_store_ps((float*)tmp, mL1);
	sum[1] = tmp[0];

	_mm_store_ps((float*)tmp, mT1);
	sum[2] = tmp[0];
}
#endif

/*
 * Reference SIMD implementation for k-fold binning sum/asum/nrm2
 */
#if defined( SCASUMI2 )
void scasumI2(int n, float complex* v, int incv, int fold, float complex* sum) {
#endif
#if defined( SCNRM2I2 )
void scnrm2I2(int n, float complex* v, int incv, float scale,
	int fold, float complex* sum) {
#endif
#if defined( CSUMI2 )
void csumI2(int n, float complex* v, int incv, int fold, float complex* sum) {
#endif

	SET_DAZ_FLAG

#ifdef MANUAL_UNROLL
	if (fold == 3) {
#		if defined( SCASUMI2 )
		scasumI2_k3(n, v, incv, sum);
#		endif
#		if defined( SCNRM2I2 )
		scnrm2I2_k3(n, v, incv, scale, sum);
#		endif
#		if defined( CSUMI2 )
		csumI2_k3(n, v, incv, sum);
#		endif
		RESET_DAZ_FLAG
		return;
	}
#endif
	
	int i, j;
	__m128 mBUF[MAX_FOLD];
	__m128 mv0, mv1;
	__m128 mq0, mq1;
	__m128 mS0, mR;
	__m128 mBLP;
	float tmp[4] __attribute__ ((aligned(16)));
	float* ptr;

	// ABSOLUTE MASK
#	if defined( SCASUMI2 )
	__m128 mAbsMask;
	SIMD_ABS_MASKS(mAbsMask);
#	endif
#	ifdef DNRM2I2
	__m128 mScale;
	mScale = _mm_set1_ps(scale);
#	endif

	// SET LAST BIT MASK: mBLP = 0.000...1 * 2^0
	SIMD_BLP_MASKS(mBLP);

	// EXPAND INITIAL SUM TO BUFFER
	for (j = 0; j < fold; j++) {
		mBUF[j] = _mm_loadu_ps((float*)(sum + j));
	}

	i = 0;


#if (CSUMI2_UNROLL_STEP >= 4)
	for (; i < n-1; i+=2, v += 2*incv) {
		mv0 = _mm_loadu_ps((float*)v);
		mv1 = _mm_loadu_ps((float*)(v + incv));

#		if defined( SCASUMI2 )
		mv0 = _mm_and_ps(mv0, mAbsMask);
		mv1 = _mm_and_ps(mv1, mAbsMask);
#		endif
#		ifdef DNRM2I
		mv0 = _mm_mul_sd(mv0, mScale);
		mv0 = _mm_mul_sd(mv0, mv0);
		mv1 = _mm_mul_ps(mv1, mScale);
		mv1 = _mm_mul_ps(mv1, mv1);
#		endif

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

#if (CSUMI2_UNROLL_STEP >= 2)
	for (; i < n; i++, v+=incv) {
		mv0 = _mm_loadu_ps((float*)v);
#		if defined( SCASUMI2 )
		mv0 = _mm_and_ps(mv0, mAbsMask);
#		endif
#		ifdef DNRM2I
		mv0 = _mm_mul_sd(mv0, mScale);
		mv0 = _mm_mul_sd(mv0, mv0);
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
#endif

	for (j = 0; j < fold; j++) {
		_mm_store_ps(tmp, mBUF[j]);
		ptr = (float*) (sum + j);
		ptr[0] = tmp[0];
		ptr[1] = tmp[1];
	}
	RESET_DAZ_FLAG
}

