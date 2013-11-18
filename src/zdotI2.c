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
#include "types.h"

#define MANUAL_UNROLL

#define UNROLL_STEP_k3 2

/*
 * Reference SIMD implementation for 3-fold binning sum/asum/nrm2/dot
 */

#ifdef MANUAL_UNROLL

#ifdef ZDOTCI2
void zdotcI2_k3
#elif defined (ZDOTUI2)
void zdotuI2_k3
#endif
(int n, double complex* v, int incv, double complex* y, int incy, double complex* sum)
{

	int i;

	__m128d mR0, mI0;
	__m128d mR2, mI2;
	__m128d mR0p, mI0p;
	__m128d mR1, mI1;
	__m128d mv0, mv1;
	__m128d my0, my1;
	__m128d mvr0, mvr1;
	double real, imag;

	// GENERATING MASK
	// mBLP = 0.000...1 * 2^0
	__m128d mBLP; // BIT IN THE LAST PLACE MASK
	SIMD_BLP_MASKD(mBLP);
	double complex tmp[1] __attribute__((aligned(16)));


	mR0 = _mm_set1_pd(ZREAL_(sum[0]));
	mI0 = _mm_set1_pd(ZIMAG_(sum[0]));

	mR1 = _mm_set1_pd(ZREAL_(sum[1]));
	mI1 = _mm_set1_pd(ZIMAG_(sum[1]));

	mR2 = _mm_set1_pd(ZREAL_(sum[2]));
	mI2 = _mm_set1_pd(ZIMAG_(sum[2]));

	i = 0;
	
#if (UNROLL_STEP_k3 >= 4)
	for (; i < n-1; i+=2, v+=2*incv, y+=2*incy) {
		mv0 = _mm_loadu_pd((double*)v);
		mv1 = _mm_loadu_pd((double*)(v + incv));

		my0 = _mm_loadu_pd((double*)y);
		my1 = _mm_loadu_pd((double*)(y + incy));

		mvr0 = _mm_shuffle_pd(mv0, mv0, _MM_SHUFFLE2(0,1));
		mvr1 = _mm_shuffle_pd(mv1, mv1, _MM_SHUFFLE2(0,1));

		mv0 = _mm_mul_pd(mv0, my0);
		mvr0 = _mm_mul_pd(mvr0, my0);
		mv1 = _mm_mul_pd(mv1, my1);
		mvr1 = _mm_mul_pd(mvr1, my1);

		// CHANGE LAST BIT TO 1
		mR0p = _mm_or_pd(mv0, mBLP);
		mI0p = _mm_or_pd(mvr0, mBLP);
		my0  = _mm_or_pd(mv1, mBLP);
		my1  = _mm_or_pd(mvr1, mBLP);

		mR0p  = _mm_add_pd(mR0p, mR0);
		mI0p  = _mm_add_pd(mI0p, mI0);
		my0   = _mm_add_pd(my0, mR0p);
		my1   = _mm_add_pd(my1, mI0p);

		mR0  = _mm_sub_pd(mR0, mR0p);
		mI0  = _mm_sub_pd(mI0, mI0p);
		mR0p = _mm_sub_pd(mR0p, my0);
		mI0p = _mm_sub_pd(mI0p, my1);

		mv0  = _mm_add_pd(mv0, mR0);
		mvr0 = _mm_add_pd(mvr0, mI0);
		mv1  = _mm_add_pd(mv1, mR0p);
		mvr1 = _mm_add_pd(mvr1, mI0p);

		mR0  = my0;
		mI0  = my1;

		//----
		mR0p = _mm_or_pd(mv0, mBLP);
		mI0p = _mm_or_pd(mvr0, mBLP);
		my0  = _mm_or_pd(mv1, mBLP);
		my1  = _mm_or_pd(mvr1, mBLP);

		mR0p  = _mm_add_pd(mR0p, mR1);
		mI0p  = _mm_add_pd(mI0p, mI1);
		my0   = _mm_add_pd(my0, mR0p);
		my1   = _mm_add_pd(my1, mI0p);

		mR1  = _mm_sub_pd(mR1, mR0p);
		mI1  = _mm_sub_pd(mI1, mI0p);
		mR0p = _mm_sub_pd(mR0p, my0);
		mI0p = _mm_sub_pd(mI0p, my1);

		mv0  = _mm_add_pd(mv0, mR1);
		mvr0 = _mm_add_pd(mvr0, mI1);
		mv1  = _mm_add_pd(mv1, mR0p);
		mvr1 = _mm_add_pd(mvr1, mI0p);

		mv0  = _mm_or_pd(mv0, mBLP);
		mv1  = _mm_or_pd(mv1, mBLP);
		mvr0 = _mm_or_pd(mvr0, mBLP);
		mvr1 = _mm_or_pd(mvr1, mBLP);

		mR1  = my0;
		mI1  = my1;


		//-----
		mR2  = _mm_add_pd(mR2, mv0);
		mI2  = _mm_add_pd(mI2, mvr0);

		mR2  = _mm_add_pd(mR2, mv1);
		mI2  = _mm_add_pd(mI2, mvr1);
	}

#endif

#if (UNROLL_STEP_k3 >= 2)
#if (UNROLL_STEP_k3 >= 4)
	if (i < n) {
		i += 1;
#else
	for (; i < n; i+=1, v += incv, y += incy) {
#endif
		mv0 = _mm_loadu_pd((double*)v);
		my0 = _mm_loadu_pd((double*)y);

		mv1 = _mm_shuffle_pd(mv0, mv0, _MM_SHUFFLE2(0,1));

		mv0 = _mm_mul_pd(mv0, my0);
		mv1 = _mm_mul_pd(mv1, my0);

		mR0p = _mm_or_pd(mv0, mBLP);
		mI0p = _mm_or_pd(mv1, mBLP);

		//-----
		mR0p  = _mm_add_pd(mR0p, mR0);
		mI0p  = _mm_add_pd(mI0p, mI0);

		mR0  = _mm_sub_pd(mR0, mR0p);
		mI0  = _mm_sub_pd(mI0, mI0p);

		mv0  = _mm_add_pd(mv0, mR0);
		mv1  = _mm_add_pd(mv1, mI0);

		mR0  = mR0p;
		mI0  = mI0p;

		my0  = _mm_or_pd(mv0, mBLP);
		my1  = _mm_or_pd(mv1, mBLP);

		my0  = _mm_add_pd(my0, mR1);
		my1  = _mm_add_pd(my1, mI1);

		mR1  = _mm_sub_pd(mR1, my0);
		mI1  = _mm_sub_pd(mI1, my1);

		mv0  = _mm_add_pd(mv0, mR1);
		mv1  = _mm_add_pd(mv1, mI1);

		mv0 = _mm_or_pd(mv0, mBLP);
		mv1 = _mm_or_pd(mv1, mBLP);

		mR1  = my0;
		mI1  = my1;

		mR2  = _mm_add_pd(mR2, mv0);
		mI2  = _mm_add_pd(mI2, mv1);
	}
#endif

		_mm_store_pd((double*)tmp, mR0);

#ifdef ZDOTUI2
		real = (ZREAL_(sum[0]) - ZIMAG_(tmp[0])) + ZREAL_(tmp[0]);
#elif defined ( ZDOTCI2 )
		real = (ZIMAG_(tmp[0]) - ZREAL_(sum[0])) + ZREAL_(tmp[0]);
#endif

		_mm_store_pd((double*)tmp, mI0);
#ifdef ZDOTCI2
		imag = (ZIMAG_(sum[0]) - ZREAL_(tmp[0])) + ZIMAG_(tmp[0]);
#elif defined ( ZDOTUI2 )
		imag = (ZIMAG_(tmp[0]) - ZIMAG_(sum[0])) + ZREAL_(tmp[0]);
#endif

		ZSET_(sum[0], real, imag);

		_mm_store_pd((double*)tmp, mR1);

#ifdef ZDOTUI2
		real = (ZREAL_(sum[1]) - ZIMAG_(tmp[0])) + ZREAL_(tmp[0]);
#elif defined ( ZDOTCI2 )
		real = (ZIMAG_(tmp[0]) - ZREAL_(sum[1])) + ZREAL_(tmp[0]);
#endif

		_mm_store_pd((double*)tmp, mI1);
#ifdef ZDOTCI2
		imag = (ZIMAG_(sum[1]) - ZREAL_(tmp[0])) + ZIMAG_(tmp[0]);
#elif defined ( ZDOTUI2 )
		imag = (ZIMAG_(tmp[0]) - ZIMAG_(sum[1])) + ZREAL_(tmp[0]);
#endif

		ZSET_(sum[1], real, imag);

		_mm_store_pd((double*)tmp, mR2);

#ifdef ZDOTUI2
		real = (ZREAL_(sum[2]) - ZIMAG_(tmp[0])) + ZREAL_(tmp[0]);
#elif defined ( ZDOTCI2 )
		real = (ZIMAG_(tmp[0]) - ZREAL_(sum[2])) + ZREAL_(tmp[0]);
#endif

		_mm_store_pd((double*)tmp, mI2);
#ifdef ZDOTCI2
		imag = (ZIMAG_(sum[2]) - ZREAL_(tmp[0])) + ZIMAG_(tmp[0]);
#elif defined ( ZDOTUI2 )
		imag = (ZIMAG_(tmp[0]) - ZIMAG_(sum[2])) + ZREAL_(tmp[0]);
#endif

		ZSET_(sum[2], real, imag);
}
#endif

/*
 * Reference SIMD implementation for k-fold binning sum/asum/nrm2
 */
#ifdef ZDOTCI2
void zdotcI2
#elif defined (ZDOTUI2)
void zdotuI2
#endif
(int n, double complex* v, int incv, double complex* y, int incy, int fold, double complex* sum)
{

	unsigned int ocsr = _mm_getcsr();
	unsigned int ncsr = ocsr | 0x8040;
	if (ncsr != ocsr) _mm_setcsr(ncsr);

#ifdef MANUAL_UNROLL
	if (fold == 3) {
#		if defined( ZDOTCI2 )
		zdotcI2_k3(n, v, incv, y, incy, sum);
#		endif
#		if defined( ZDOTUI2 )
		zdotuI2_k3(n, v, incv, y, incy, sum);
#		endif
		if (ncsr != ocsr) _mm_setcsr(ocsr);
		return;
	}
#endif
	
	int i, j;
	__m128d mR[MAX_FOLD];
	__m128d mI[MAX_FOLD];
	__m128d mv0, mv1;
	__m128d my0;
	__m128d mq0, mq1;
	__m128d mR0, mI0;
	double complex tmp[1] __attribute__((aligned(16)));
	double real, imag;

	__m128d mBLP;

	// ABSOLUTE MASK

	// SET LAST BIT MASK: mBLP = 0.000...1 * 2^0
	SIMD_BLP_MASKD(mBLP);

	// EXPAND INITIAL SUM TO BUFFER
	for (j = 0; j < fold; j++) {
		mR[j] = _mm_set1_pd(ZREAL_(sum[j]));
		mI[j] = _mm_set1_pd(ZIMAG_(sum[j]));
	}

	i = 0;


	for (; i < n-1; i+=2, v += 2*incv) {
		mv0 = _mm_loadu_pd((double*)v);
		my0 = _mm_loadu_pd((double*)y);

		mv1 = _mm_shuffle_pd(mv0, mv0, _MM_SHUFFLE2(0,1));

		mv0 = _mm_mul_pd(mv0, my0);
		mv1 = _mm_mul_pd(mv1, my0);


		for (j = 0; j < fold - 1; j++) {
			mq0     = _mm_or_pd(mv0, mBLP);
			mq1     = _mm_or_pd(mv1, mBLP);

			mR0     = mR[j];
			mI0     = mI[j];

			mR[j] = _mm_add_pd(mR0, mq0);
			mI[j] = _mm_add_pd(mI0, mq1);

			mR0     = _mm_sub_pd(mR0, mR[j]);
			mI0     = _mm_sub_pd(mI0, mI[j]);

			mv0     = _mm_add_pd(mv0, mR0);
			mv1     = _mm_add_pd(mv1, mI0);
		}

		mv0   = _mm_or_pd(mv0, mBLP);
		mv1   = _mm_or_pd(mv1, mBLP);
		mR[j] = _mm_add_pd(mR[j], mv0);
		mI[j] = _mm_add_pd(mI[j], mv1);
	}
	
	for (j = 0; j < fold; j++) {
		_mm_store_pd((double*)tmp, mR[j]);

#ifdef ZDOTUI2
		real = (ZREAL_(sum[j]) - ZIMAG_(tmp[0])) + ZREAL_(tmp[0]);
#elif defined ( ZDOTCI2 )
		real = (ZIMAG_(tmp[0]) - ZREAL_(sum[j])) + ZREAL_(tmp[0]);
#endif

		_mm_store_pd((double*)tmp, mI[j]);
#ifdef ZDOTCI2
		imag = (ZIMAG_(sum[i]) - ZREAL_(tmp[0])) + ZIMAG_(tmp[0]);
#elif defined ( ZDOTUI2 )
		imag = (ZIMAG_(tmp[0]) - ZIMAG_(sum[j])) + ZREAL_(tmp[0]);
#endif

		ZSET_(sum[j], real, imag);
	}
	if (ncsr != ocsr) _mm_setcsr(ocsr);
}

