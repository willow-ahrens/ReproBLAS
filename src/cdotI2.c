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
#include "types.h"

#define MANUAL_UNROLL

#define UNROLL_STEP_NR_k3 2

/*
 * Reference SIMD implementation for 3-fold binning sum/asum/nrm2/dot
 */

#ifdef MANUAL_UNROLL

#ifdef CDOTCI2
void cdotcI2_k3
#elif defined (CDOTUI2)
void cdotuI2_k3
#endif
(int n, float complex* v, int incv, float complex* y, int incy, float complex* sum)
{

	int i;

	__m128 mR0, mI0;
	__m128 mR2, mI2;
	__m128 mR0p, mI0p;
	__m128 mR1, mI1;
	__m128 mR1p, mI1p;
	__m128 mv0, mv1;
	__m128 my0, my1;
	__m128 mvr0, mvr1;
	float real, imag;

	// GENERATING MASK
	// mBLP = 0.000...1 * 2^0
	__m128 mBLP; // BIT IN THE LAST PLACE MASK
	SIMD_BLP_MASKS(mBLP);
	float complex sret[2] __attribute__((aligned(16)));

	mR0 = _mm_set1_ps(CREAL_(sum[0]));
	mI0 = _mm_set1_ps(CIMAG_(sum[0]));

	mR0p = mR0;
	mI0p = mI0;

	mR1 = _mm_set1_ps(CREAL_(sum[1]));
	mI1 = _mm_set1_ps(CIMAG_(sum[1]));

	mR1p = mR1;
	mI1p = mI1;

	mR2 = _mm_set1_ps(CREAL_(sum[2]));
	mI2 = _mm_set1_ps(CIMAG_(sum[2]));

	i = 0;
	
#if (UNROLL_STEP_NR_k3 >= 4)
	for (; i < n-1; i+=2, v+=2*incv, y += 2 * incy) {
		mv0 = _mm_loadu_ps((float*)v);
		mv1 = _mm_loadu_ps((float*)(v + incv));

		my0 = _mm_loadu_ps((float*)y);
		my1 = _mm_loadu_ps((float*)(y + incv));

		mvr0 = _mm_shuffle_ps(mv0, mv0, _MM_SHUFFLE(0,1,0,1));
		mvr1 = _mm_shuffle_ps(mv1, mv1, _MM_SHUFFLE(0,1,0,1));

		mv0 = _mm_mul_ps(mv0, my0);
		mvr0 = _mm_mul_ps(mvr0, my0);
		mv1 = _mm_mul_ps(mv1, my1);
		mvr1 = _mm_mul_ps(mvr1, my1);

		// CHANGE LAST BIT TO 1
		mR0p = _mm_or_ps(mv0, mBLP);
		mI0p = _mm_or_ps(mv1, mBLP);

		//-----
		mR0p  = _mm_add_ps(mR0p, mR0);
		mR0  = _mm_sub_ps(mR0, mR0p);
		mv0  = _mm_add_ps(mv0, mR0);
		mR0  = mR0p;

		mR1p = _mm_or_ps(mv0, mBLP);
		mR1p  = _mm_add_ps(mR1p, mR1);
		mR1  = _mm_sub_ps(mR1, mR1p);
		mv0  = _mm_add_ps(mv0, mR1);
		mv0 = _mm_or_ps(mv0, mBLP);
		mR1  = mR1p;

		mR2  = _mm_add_ps(mR2, mv0);
		//----

		mI0p  = _mm_add_ps(mI0p, mI0);
		mI0  = _mm_sub_ps(mI0, mI0p);
		mv1  = _mm_add_ps(mv1, mI0);
		mI0  = mI0p;

		mI1p = _mm_or_ps(mv1, mBLP);
		mI1p  = _mm_add_ps(mI1p, mI1);
		mI1  = _mm_sub_ps(mI1, mI1p);
		mv1  = _mm_add_ps(mv1, mI1);
		mv1 = _mm_or_ps(mv1, mBLP);
		mI1  = mI1p;

		mI2  = _mm_add_ps(mI2, mv1);
		//-----
	}

	mR0p = _mm_loadu_ps((float*)sum);
	mR1p = _mm_loadu_ps((float*)(sum+1));
	mI0p = _mm_loadu_ps((float*)(sum+2));

	mI0 = _mm_sub_ps(mI0, mR0p);
	mI1 = _mm_sub_ps(mI1, mR1p);
	mI2 = _mm_sub_ps(mI2, mI0p);

	mR0 = _mm_add_ps(mR0, mI0);
	mR1 = _mm_add_ps(mR1, mI1);
	mR2 = _mm_add_ps(mR2, mI2);

	mR0p = mR0;
	mR1p = mR1;
#endif

#if (UNROLL_STEP_NR_k3 >= 2)
#if (UNROLL_STEP_NR_k3 >= 4)
	if (i < n) {
		i += 1;
#else
	for (; i < n; i+=1) {
#endif
		mv0 = _mm_loadu_ps((float*)v);
		my0 = _mm_loadu_ps((float*)y);

		mv1 = _mm_shuffle_ps(mv0, mv0, _MM_SHUFFLE(0,1,0,1));

		mv0 = _mm_mul_ps(mv0, my0);
		mv1 = _mm_mul_ps(mv1, my0);

		mR0p = _mm_or_ps(mv0, mBLP);
		mI0p = _mm_or_ps(mv1, mBLP);

		//-----
		mR0p  = _mm_add_ps(mR0p, mR0);
		mI0p  = _mm_add_ps(mI0p, mI0);

		mR0  = _mm_sub_ps(mR0, mR0p);
		mI0  = _mm_sub_ps(mI0, mI0p);

		mv0  = _mm_add_ps(mv0, mR0);
		mv1  = _mm_add_ps(mv1, mI0);

		mR0  = mR0p;
		mI0  = mI0p;

		mR1p = _mm_or_ps(mv0, mBLP);
		mI1p = _mm_or_ps(mv1, mBLP);

		mR1p  = _mm_add_ps(mR1p, mR1);
		mI1p  = _mm_add_ps(mI1p, mI1);

		mR1  = _mm_sub_ps(mR1, mR1p);
		mI1  = _mm_sub_ps(mI1, mI1p);

		mv0  = _mm_add_ps(mv0, mR1);
		mv1  = _mm_add_ps(mv1, mI1);

		mv0 = _mm_or_ps(mv0, mBLP);
		mv1 = _mm_or_ps(mv1, mBLP);

		mR1  = mR1p;
		mI1  = mI1p;

		mR2  = _mm_add_ps(mR2, mv0);
		mI2  = _mm_add_ps(mI2, mv1);
		//----
		v += incv;
		y += incy;
	}
#endif

	_mm_store_ps((float*) sret, mR0);
#ifdef CDOTUI2
	real = (CREAL_(sum[0]) - CIMAG_(sret[0])) + CREAL_(sret[0]);
#elif defined ( CDOTCI2 )
	real = (CREAL_(sret[0])- CREAL_(sum[0])) + CIMAG_(sret[0]);
#endif

	_mm_store_ps((float*) sret, mI0);
#ifdef CDOTCI2
	imag = (CIMAG_(sum[0]) - CREAL_(sret[0])) + CIMAG_(sret[0]);
#elif defined ( CDOTUI2 )
	imag = (CIMAG_(sret[0]) - CIMAG_(sum[0])) + CREAL_(sret[0]);
#endif
	CSET_(sum[0], real, imag);
	//----

	_mm_store_ps((float*) sret, mR1);
#ifdef CDOTUI2
	real = (CREAL_(sum[1]) - CIMAG_(sret[0])) + CREAL_(sret[0]);
#elif defined ( CDOTCI2 )
	real = (CREAL_(sret[0])- CREAL_(sum[1])) + CIMAG_(sret[0]);
#endif

	_mm_store_ps((float*) sret, mI1);
#ifdef CDOTCI2
	imag = (CIMAG_(sum[1]) - CREAL_(sret[0])) + CIMAG_(sret[0]);
#elif defined ( CDOTUI2 )
	imag = (CIMAG_(sret[0]) - CIMAG_(sum[1])) + CREAL_(sret[0]);
#endif
	CSET_(sum[1], real, imag);
	//----

	_mm_store_ps((float*) sret, mR2);
#ifdef CDOTUI2
	real = (CREAL_(sum[2]) - CIMAG_(sret[0])) + CREAL_(sret[0]);
#elif defined ( CDOTCI2 )
	real = (CREAL_(sret[0])- CREAL_(sum[2])) + CIMAG_(sret[0]);
#endif

	_mm_store_ps((float*) sret, mI2);
#ifdef CDOTCI2
	imag = (CIMAG_(sum[2]) - CREAL_(sret[0])) + CIMAG_(sret[0]);
#elif defined ( CDOTUI2 )
	imag = (CIMAG_(sret[0]) - CIMAG_(sum[2])) + CREAL_(sret[0]);
#endif
	CSET_(sum[2], real, imag);
	//----
}
#endif

/*
 * Reference SIMD implementation for k-fold binning sum/asum/nrm2
 */
#ifdef CDOTCI2
void cdotcI2
#elif defined (CDOTUI2)
void cdotuI2
#endif
(int n, float complex* v, int incv, float complex* y, int incy, int fold, float complex* sum)
{

	unsigned int ocsr = _mm_getcsr();
	unsigned int ncsr = ocsr | 0x8040;
	if (ncsr != ocsr) _mm_setcsr(ncsr);

#ifdef MANUAL_UNROLL
	if (fold == 3) {
#		if defined( CDOTCI2 )
		cdotcI2_k3(n, v, incv, y, incy, sum);
#		endif
#		if defined( CDOTUI2 )
		cdotuI2_k3(n, v, incv, y, incy, sum);
#		endif
		if (ncsr != ocsr) _mm_setcsr(ocsr);
		return;
	}
#endif
	
	int i, j;
	__m128 mR[MAX_FOLD];
	__m128 mI[MAX_FOLD];
	__m128 mv0, mv1;
	__m128 my0;
	__m128 mq0, mq1;
	__m128 mR0, mI0;
	float complex sret[2] __attribute__((aligned(16)));
	float real, imag;

	__m128 mBLP;

	// ABSOLUTE MASK

	// SET LAST BIT MASK: mBLP = 0.000...1 * 2^0
	SIMD_BLP_MASKS(mBLP);

	// EXPAND INITIAL SUM TO BUFFER
	for (j = 0; j < fold; j++) {
		mR[j] = _mm_set1_ps(CREAL_(sum[j]));
		mI[j] = _mm_set1_ps(CIMAG_(sum[j]));
	}

	i = 0;


	for (; i < n-1; i+=2, v += 2*incv) {
		mv0 = _mm_loadu_ps((float*)v);
		my0 = _mm_loadu_ps((float*)y);

		mv1 = _mm_shuffle_ps(mv0, mv0, _MM_SHUFFLE(0,1,0,1));

		mv0 = _mm_mul_ps(mv0, my0);
		mv1 = _mm_mul_ps(mv1, my0);


		for (j = 0; j < fold - 1; j++) {
			mq0     = _mm_or_ps(mv0, mBLP);
			mq1     = _mm_or_ps(mv1, mBLP);

			mR0     = mR[j];
			mI0     = mI[j];

			mR[j] = _mm_add_ps(mR0, mq0);
			mI[j] = _mm_add_ps(mI0, mq1);

			mR0     = _mm_sub_ps(mR0, mR[j]);
			mI0     = _mm_sub_ps(mI0, mI[j]);

			mv0     = _mm_add_ps(mv0, mR0);
			mv1     = _mm_add_ps(mv1, mI0);
		}

		mv0   = _mm_or_ps(mv0, mBLP);
		mv1   = _mm_or_ps(mv1, mBLP);
		mR[j] = _mm_add_ps(mR[j], mv0);
		mI[j] = _mm_add_ps(mI[j], mv1);
	}
	
	for (j = 0; j < fold; j++) {
		_mm_store_ps((float*) sret, mR[j]);
#ifdef CDOTUI2
		real = (CREAL_(sum[j]) - CREAL_(sret[0])) + CIMAG_(sret[0]);
#elif defined ( CDOTCI2 )
		real = (CREAL_(sret[0])- CREAL_(sum[j])) + CIMAG_(sret[0]);
#endif

		_mm_store_ps((float*) sret, mI[j]);
#ifdef CDOTCI2
		imag = (CIMAG_(sum[j]) - CIMAG_(sret[0])) + CREAL_(sret[0]);
#elif defined ( CDOTUI2 )
		imag = (CIMAG_(sret[0]) - CIMAG_(sum[j])) + CREAL_(sret[0]);
#endif
		CSET_(sum[j], real, imag);
	}
	if (ncsr != ocsr) _mm_setcsr(ocsr);
}

