/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <emmintrin.h>

#include "config.h"
#include "Common/Common.h"
#include "zIBlas1.h"
#include "rblas1.h"

#define MANUAL_UNROLL

#define UNROLL_STEP_NR_k3 2

/*
 * Reference SIMD implementation for 3-fold binning sum/asum/nrm2/dot
 */

#ifdef MANUAL_UNROLL

#if defined( DZASUMI2 )
void dzasumI2_k3(int n, double complex* v, int incv, double complex* sum) {
#endif
#if defined( ZSUMI2 )
void zsumI2_k3(int n, double complex* v, int incv, double complex* sum) {
#endif
#if defined( DZNRM2I2 )
void dznrm2I2_k3(int n, double complex* v, int incv, double scale, double complex* sum) {
#endif

	int i;

	__m128d mS1, mS2;
	__m128d mT1, mT2;
	__m128d mS1p, mS2p;
	__m128d mL1, mL2;
	__m128d mL1p, mL2p;
	__m128d mv0, mv1;

	// GENERATING MASK
	// mBLP = 0.000...1 * 2^0
	__m128d mBLP; // BIT IN THE LAST PLACE MASK
	SSE_BLP_MASKD(mBLP);

#if defined( DZASUMI2 )
	__m128d mAbsMask;
	SSE_ABS_MASKD(mAbsMask);
#endif
#ifdef DZNRM2I2
	__m128d mScale;
	mScale = _mm_set1_pd(scale);
#endif

	mS1 = _mm_loadu_pd((double*)sum);
	mS2 = mS1;

	mS1p = mS1;
	mS2p = mS2;

	mT1 = _mm_loadu_pd((double*)(sum + 2));
	mT2 = mT1;

	mL1 = _mm_loadu_pd((double*)(sum + 1));
	mL2 = mL1;

	mL1p = mL1;
	mL2p = mL2;

	i = 0;
	
#if (UNROLL_STEP_NR_k3 >= 4)
	for (; i < n-1; i+=2, v+=2*incv) {
		mv0 = _mm_loadu_pd((double*)v);
		mv1 = _mm_loadu_pd((double*)(v + incv));

#		ifdef DZASUMI2
		// ABSOLUTE VALUE
		mv0 = _mm_and_pd(mv0, mAbsMask);
		mv1 = _mm_and_pd(mv1, mAbsMask);
#		endif
#		ifdef DZNRM2I2
		mv0 = _mm_mul_pd(mv0, mScale);
		mv0 = _mm_mul_pd(mv0, mv0);
		mv1 = _mm_mul_pd(mv1, mScale);
		mv1 = _mm_mul_pd(mv1, mv1);
#endif
		// CHANGE LAST BIT TO 1
		mS1p = _mm_or_pd(mv0, mBLP);
		mS2p = _mm_or_pd(mv1, mBLP);

		//-----
		mS1p  = _mm_add_pd(mS1p, mS1);
		mS1  = _mm_sub_pd(mS1, mS1p);
		mv0  = _mm_add_pd(mv0, mS1);
		mS1  = mS1p;

		mL1p = _mm_or_pd(mv0, mBLP);
		mL1p  = _mm_add_pd(mL1p, mL1);
		mL1  = _mm_sub_pd(mL1, mL1p);
		mv0  = _mm_add_pd(mv0, mL1);
		mv0 = _mm_or_pd(mv0, mBLP);
		mL1  = mL1p;

		mT1  = _mm_add_pd(mT1, mv0);
		//----

		mS2p  = _mm_add_pd(mS2p, mS2);
		mS2  = _mm_sub_pd(mS2, mS2p);
		mv1  = _mm_add_pd(mv1, mS2);
		mS2  = mS2p;

		mL2p = _mm_or_pd(mv1, mBLP);
		mL2p  = _mm_add_pd(mL2p, mL2);
		mL2  = _mm_sub_pd(mL2, mL2p);
		mv1  = _mm_add_pd(mv1, mL2);
		mv1 = _mm_or_pd(mv1, mBLP);
		mL2  = mL2p;

		mT2  = _mm_add_pd(mT2, mv1);
		//-----
	}

	mS1p = _mm_loadu_pd((double*)sum);
	mL1p = _mm_loadu_pd((double*)(sum+1));
	mS2p = _mm_loadu_pd((double*)(sum+2));

	mS2 = _mm_sub_pd(mS2, mS1p);
	mL2 = _mm_sub_pd(mL2, mL1p);
	mT2 = _mm_sub_pd(mT2, mS2p);

	mS1 = _mm_add_pd(mS1, mS2);
	mL1 = _mm_add_pd(mL1, mL2);
	mT1 = _mm_add_pd(mT1, mT2);

	mS1p = mS1;
	mL1p = mL1;
#endif

#if (UNROLL_STEP_NR_k3 >= 2)
#if (UNROLL_STEP_NR_k3 >= 4)
	if (i < n) {
		i += 1;
#else
	for (; i < n; i+=1) {
#endif
		mv0 = _mm_loadu_pd((double*)v);

#if defined( DZASUMI2 )
		mv0 = _mm_and_pd(mv0, mAbsMask);
#endif
#		ifdef DZNRM2I2
		mv0 = _mm_mul_pd(mv0, mScale);
		mv0 = _mm_mul_pd(mv0, mv0);
#endif
		mS1p = _mm_or_pd(mv0, mBLP);

		//-----
		mS1p  = _mm_add_pd(mS1p, mS1);
		mS1  = _mm_sub_pd(mS1, mS1p);
		mv0  = _mm_add_pd(mv0, mS1);
		mS1  = mS1p;

		mL1p = _mm_or_pd(mv0, mBLP);

		mL1p  = _mm_add_pd(mL1p, mL1);
		mL1  = _mm_sub_pd(mL1, mL1p);
		mv0  = _mm_add_pd(mv0, mL1);
		mv0 = _mm_or_pd(mv0, mBLP);
		mL1  = mL1p;

		mT1  = _mm_add_pd(mT1, mv0);
		//----
		v += incv;
	}
#endif

	_mm_storeu_pd((double*)sum, mS1);
	_mm_storeu_pd((double*)(sum + 1), mL1);
	_mm_storeu_pd((double*)(sum + 2), mT1);
}
#endif

/*
 * Reference SIMD implementation for k-fold binning sum/asum/nrm2
 */
#if defined( DZASUMI2 )
void dzasumI2(int n, double complex* v, int incv, int fold, double complex* sum) {
#endif
#if defined( DZNRM2I2 )
void dznrm2I2(int n, double complex* v, int incv, double scale,
	int fold, double complex* sum) {
#endif
#if defined( ZSUMI2 )
void zsumI2(int n, double complex* v, int incv, int fold, double complex* sum) {
#endif

	SET_DAZ_FLAG

#ifdef MANUAL_UNROLL
	if (fold == 3) {
#		if defined( DZASUMI2 )
		dzasumI2_k3(n, v, incv, sum);
#		endif
#		if defined( DZNRM2I2 )
		dznrm2I2_k3(n, v, incv, scale, sum);
#		endif
#		if defined( ZSUMI2 )
		zsumI2_k3(n, v, incv, sum);
#		endif
		RESET_DAZ_FLAG
		return;
	}
#endif
	
	int i, j;
	__m128d mBUF[MAX_FOLD];
	__m128d mv0, mv1;
	__m128d mq0, mq1;
	__m128d mS0, mR;
	__m128d mBLP;

	// ABSOLUTE MASK
#	if defined( DZASUMI2 )
	__m128d mAbsMask;
	SSE_ABS_MASKD(mAbsMask);
#	endif
#	ifdef DNRM2I2
	__m128d mScale;
	mScale = _mm_set1_pd(scale);
#	endif

	// SET LAST BIT MASK: mBLP = 0.000...1 * 2^0
	SSE_BLP_MASKD(mBLP);

	// EXPAND INITIAL SUM TO BUFFER
	for (j = 0; j < fold; j++) {
		mBUF[j] = _mm_loadu_pd((double*)(sum + j));
	}

	i = 0;


#if (DSUMI2_UNROLL_STEP >= 4)
	for (; i < n-1; i+=2, v += 2*incv) {
		mv0 = _mm_loadu_pd((double*)v);
		mv1 = _mm_loadu_pd((double*)(v + incv));

#		if defined( DZASUMI2 )
		mv0 = _mm_and_pd(mv0, mAbsMask);
		mv1 = _mm_and_pd(mv1, mAbsMask);
#		endif
#		ifdef DNRM2I
		mv0 = _mm_mul_sd(mv0, mScale);
		mv0 = _mm_mul_sd(mv0, mv0);
		mv1 = _mm_mul_pd(mv1, mScale);
		mv1 = _mm_mul_pd(mv1, mv1);
#		endif

		for (j = 0; j < fold - 1; j++) {
			mq0     = _mm_or_pd(mv0, mBLP);
			mq1     = _mm_or_pd(mv1, mBLP);

			mS0     = mBUF[j];
			mR     = _mm_add_pd(mS0, mq0);
			mBUF[j] = _mm_add_pd(mR, mq1);

			mS0     = _mm_sub_pd(mS0, mR);
			mR     = _mm_sub_pd(mR, mBUF[j]);

			mv0     = _mm_add_pd(mv0, mS0);
			mv1     = _mm_add_pd(mv1, mR);
		}

		mv0     = _mm_or_pd(mv0, mBLP);
		mv1     = _mm_or_pd(mv1, mBLP);
		mBUF[j] = _mm_add_pd(mBUF[j], mv0);
		mBUF[j] = _mm_add_pd(mBUF[j], mv1);
	}
	
#endif

#if (DSUMI2_UNROLL_STEP >= 2)
	for (; i < n; i++, v+=incv) {
		mv0 = _mm_loadu_pd((double*)v);
#		if defined( DZASUMI2 )
		mv0 = _mm_and_pd(mv0, mAbsMask);
#		endif
#		ifdef DNRM2I
		mv0 = _mm_mul_sd(mv0, mScale);
		mv0 = _mm_mul_sd(mv0, mv0);
#		endif

		for (j = 0; j < fold - 1; j++) {
			mq1     = _mm_or_pd(mv0, mBLP);
			mR     = mBUF[j];
			mBUF[j] = _mm_add_pd(mBUF[j], mq1);
			mR     = _mm_sub_pd(mR, mBUF[j]);
			mv0     = _mm_add_pd(mv0, mR);
		}

		mv0     = _mm_or_pd(mv0, mBLP);
		mBUF[j] = _mm_add_pd(mBUF[j], mv0);
	}
#endif

	for (j = 0; j < fold; j++) {
		_mm_storeu_pd((double*)(sum + j), mBUF[j]);
	}

	RESET_DAZ_FLAG
}

