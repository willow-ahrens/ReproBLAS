/*
 *  Created   13/10/25   H.D. Nguyen
 */

// REDIRECT TO SISD IMPLEMENTATION IF NO SSE SUPPORT
#ifndef __SSE2__
#include "dIBlas1_SISD.c"
#else

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <emmintrin.h>

#include "config.h"
#include "Common/Common.h"
#include "dIBlas1.h"

/*
 * SUMMATION IN INDEXEDFP FORMAT
 * REFERENCE IMPLEMENTATION
 */

extern void dsumI2_k3(int, double*, int, double*);

#define MANUAL_UNROLL

#define UNROLL_STEP_NR_k3 8

/*
 * Reference SIMD implementation for 3-fold binning sum/asum/nrm2/dot
 */

#ifdef MANUAL_UNROLL

#if defined( DASUMI2 )
void dasumI2_k3(int n, double* v, int incv, double* sum) {
#elif defined( DNRM2I2 )
void dnrm2I2_k3(int n, double* v, int incv, double scale, double* sum) {
#elif defined( DSUMI2 )
void dsumI2_k3(int n, double* v, int incv, double* sum) {
#elif defined( DDOTI2 )
void ddotI2_k3(int n, double* v, int incv, double* y, int incy, double* sum) {
#endif
	int i;

	__m128d mS1, mS2;
	__m128d mT1, mT2;
	__m128d mS1p, mS2p;
	__m128d mL1, mL2;
	__m128d mL1p, mL2p;
	__m128d mv0, mv1;
#if (UNROLL_STEP_NR_k3 >= 8)
	__m128d mv2, mv3;
#endif
	__m128d mBLP; // BIT IN THE LAST PLACE MASK
	double tmp[2] __attribute__((aligned(16)));

	// GENERATING MASK
	// mBLP = 0.000...1 * 2^0
	SSE_BLP_MASKD(mBLP);

	// ABSOLUTE MASK
#if defined( DASUMI2 )
	__m128d mAbsMask;
	SSE_ABS_MASKD(mAbsMask);
#elif defined( DNRM2I2 )
	__m128d mScale;
	mScale = _mm_set1_pd(scale);
#elif defined( DDOTI2 )
	__m128d my0, my1;
#if (UNROLL_STEP_NR_k3 >= 8)
	__m128d my2, my3;
#endif
#endif

	mS1 = _mm_load1_pd(sum);
	mS2 = mS1;

	mS1p = mS1;
	mS2p = mS2;

	mT1 = _mm_load1_pd(sum + 2);
	mT2 = mT1;

	mL1 = _mm_load1_pd(sum + 1);
	mL2 = mL1;

	mL1p = mL1;
	mL2p = mL2;

	i = 0;
	
#if (UNROLL_STEP_NR_k3 >= 8)
	for (; i < n-7; i+=8, v+=8) {
		mv0 = _mm_load_pd(v + 0);
		mv1 = _mm_load_pd(v + 2);
		mv2 = _mm_load_pd(v + 4);
		mv3 = _mm_load_pd(v + 6);

#if defined( DASUMI2 )
		mv0 = _mm_and_pd(mv0, mAbsMask);
		mv1 = _mm_and_pd(mv1, mAbsMask);
		mv2 = _mm_and_pd(mv2, mAbsMask);
		mv3 = _mm_and_pd(mv3, mAbsMask);
#elif defined( DNRM2I2 )
		mv0 = _mm_mul_pd(mv0, mScale);
		mv0 = _mm_mul_pd(mv0, mv0);
		mv1 = _mm_mul_pd(mv1, mScale);
		mv1 = _mm_mul_pd(mv1, mv1);
		mv2 = _mm_mul_pd(mv2, mScale);
		mv2 = _mm_mul_pd(mv2, mv2);
		mv3 = _mm_mul_pd(mv3, mScale);
		mv3 = _mm_mul_pd(mv3, mv3);
#elif defined( DDOTI2 )
		my0 = _mm_loadu_pd(y);
		my1 = _mm_loadu_pd(y + 2);
		my2 = _mm_loadu_pd(y + 4);
		my3 = _mm_loadu_pd(y + 6);
		mv0 = _mm_mul_pd(mv0, my0);
		mv1 = _mm_mul_pd(mv1, my1);
		mv2 = _mm_mul_pd(mv2, my2);
		mv3 = _mm_mul_pd(mv3, my3);
		y += 8;
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
		mL1p = _mm_add_pd(mL1p, mL1);
		mL1  = _mm_sub_pd(mL1, mL1p);
		mv0  = _mm_add_pd(mv0, mL1);
		mv0  = _mm_or_pd(mv0, mBLP);
		mL1  = mL1p;

		mT1  = _mm_add_pd(mT1, mv0);
		//----

		mS2p = _mm_add_pd(mS2p, mS2);
		mS2  = _mm_sub_pd(mS2, mS2p);
		mv1  = _mm_add_pd(mv1, mS2);
		mS2  = mS2p;

		mL2p = _mm_or_pd(mv1, mBLP);
		mL2p = _mm_add_pd(mL2p, mL2);
		mL2  = _mm_sub_pd(mL2, mL2p);
		mv1  = _mm_add_pd(mv1, mL2);
		mv1  = _mm_or_pd(mv1, mBLP);
		mL2  = mL2p;

		mT2  = _mm_add_pd(mT2, mv1);

		//======
		mS1p = _mm_or_pd(mv2, mBLP);
		mS2p = _mm_or_pd(mv3, mBLP);
		//----
		mS1p = _mm_add_pd(mS1p, mS1);
		mS1  = _mm_sub_pd(mS1, mS1p);
		mv2  = _mm_add_pd(mv2, mS1);
		mS1  = mS1p;

		mL1p = _mm_or_pd(mv2, mBLP);
		mL1p = _mm_add_pd(mL1p, mL1);
		mL1  = _mm_sub_pd(mL1, mL1p);
		mv2  = _mm_add_pd(mv2, mL1);
		mv2  = _mm_or_pd(mv2, mBLP);
		mL1  = mL1p;

		mT1  = _mm_add_pd(mT1, mv2);
		//----

		mS2p = _mm_add_pd(mS2p, mS2);
		mS2  = _mm_sub_pd(mS2, mS2p);
		mv3  = _mm_add_pd(mv3, mS2);
		mS2  = mS2p;

		mL2p = _mm_or_pd(mv3, mBLP);
		mL2p = _mm_add_pd(mL2p, mL2);
		mL2  = _mm_sub_pd(mL2, mL2p);
		mv3  = _mm_add_pd(mv3, mL2);
		mv3  = _mm_or_pd(mv3, mBLP);
		mL2  = mL2p;

		mT2  = _mm_add_pd(mT2, mv3);
		//----
	}
#endif

#if (UNROLL_STEP_NR_k3 >= 4)
	for (; i < n-3; i+=4, v+=4) {
		mv0 = _mm_load_pd(v);
		mv1 = _mm_load_pd(v + 2);

#if defined( DASUMI2 )
		mv0 = _mm_and_pd(mv0, mAbsMask);
		mv1 = _mm_and_pd(mv1, mAbsMask);
#elif defined( DNRM2I2 )
		mv0 = _mm_mul_pd(mv0, mScale);
		mv0 = _mm_mul_pd(mv0, mv0);
		mv1 = _mm_mul_pd(mv1, mScale);
		mv1 = _mm_mul_pd(mv1, mv1);
#elif defined( DDOTI2 )
		my0 = _mm_loadu_pd(y);
		my1 = _mm_loadu_pd(y + 2);
		mv0 = _mm_mul_pd(mv0, my0);
		mv1 = _mm_mul_pd(mv1, my1);
		y += 4;
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

	mS1p = _mm_set1_pd(sum[0]);
	mL1p = _mm_set1_pd(sum[1]);
	mS2p = _mm_set1_pd(sum[2]);

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
	if (i < n -1) {
		i += 2;
#else
	for (; i < n-1; i+=2) {
#endif
		mv0 = _mm_load_pd(v);
#if defined( DASUMI2 )
		mv0 = _mm_and_pd(mv0, mAbsMask);
#elif defined( DNRM2I2 )
		mv0 = _mm_mul_pd(mv0, mScale);
		mv0 = _mm_mul_pd(mv0, mv0);
#elif defined( DDOTI2 )
		my0 = _mm_loadu_pd(y);
		mv0 = _mm_mul_pd(mv0, my0);
		y += 2;
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
		v += 2;
	}
#endif

#if (UNROLL_STEP_NR_k3 >= 2)
	_mm_store_pd(tmp, mS1);
	sum[0] = tmp[1] - sum[0];
	_mm_store_pd(tmp, mL1);
	sum[1] = tmp[1] - sum[1];
	_mm_store_pd(tmp, mT1);
	sum[2] = tmp[1] - sum[2];
#endif
	// last one
#if (UNROLL_STEP_NR_k3 >= 2)
	if (i < n) {
#else
	for (; i < n; i++, v++) {
#endif
		mv0 = _mm_load_sd(v);
#		if defined( DASUMI2 )
		mv0 = _mm_and_pd(mv0, mAbsMask);
#		elif defined( DNRM2I2 )
		mv0 = _mm_mul_sd(mv0, mScale);
		mv0 = _mm_mul_sd(mv0, mv0);
#		elif defined( DDOTI2 )
		my0 = _mm_load_sd(y);
		mv0 = _mm_mul_sd(mv0, my0);
		y++;
#		endif

		mS1p  = _mm_or_pd(mv0, mBLP);

		//-----
		mS1p  = _mm_add_sd(mS1p, mS1);
		mS1   = _mm_sub_sd(mS1, mS1p);
		mv0   = _mm_add_sd(mv0, mS1);
		mS1   = mS1p;

		mL1p = _mm_or_pd(mv0, mBLP);

		mL1p  = _mm_add_sd(mL1p, mL1);
		mL1  = _mm_sub_sd(mL1, mL1p);
		mv0  = _mm_add_sd(mv0, mL1);
		mv0 = _mm_or_pd(mv0, mBLP);
		mL1  = _mm_move_sd(mL1, mL1p);

		mT1  = _mm_add_sd(mT1, mv0);
		//----
	}
	// TODO SSE3 horizontal add _mm_hadd_pd
	/*
	mS2 = _mm_load_sd(sum);
	mL2 = _mm_load_sd(sum + 1);
	mT2 = _mm_load_sd(sum + 2);
	
	mS1 = _mm_sub_pd(mS1, mS2);
	mL1 = _mm_sub_pd(mL1, mL2);
	mT1 = _mm_sub_pd(mT1, mT2);

#ifdef __SSE3___
	mS1 = _mm_hadd_pd(mS1, mS1);
	_mm_store_sd(sum, mS1);

	mL1 = _mm_hadd_pd(mL1, mL1);
	_mm_store_sd(sum + 1, mL1);

	mT1 = _mm_hadd_pd(mT1, mT1);
	_mm_store_sd(sum + 2, mT1);
#else
	_mm_store_pd(tmp, mS1);
	sum[0] = tmp[0] + tmp[1];
	_mm_store_pd(tmp, mL1);
	sum[1] = tmp[0] + tmp[1];
	_mm_store_pd(tmp, mT1);
	sum[2] = tmp[0] + tmp[1];
#endif
	*/

#if (UNROLL_STEP_NR_k3 >= 2)
	_mm_store_sd(tmp, mS1);
	sum[0] += tmp[0];
	_mm_store_sd(tmp, mL1);
	sum[1] += tmp[0];
	_mm_store_sd(tmp, mT1);
	sum[2] += tmp[0];
#else
	_mm_store_sd(sum, mS1);
	_mm_store_sd(sum+1, mL1);
	_mm_store_sd(sum+2, mT1);
#endif
}
#endif

/*
 * Reference SIMD implementation for k-fold binning sum/asum/nrm2/dot
 */
#if defined( DASUMI2 )
void dasumI2(int n, double* v, int incv, int fold, double* sum) {
#elif defined( DNRM2I2 )
void dnrm2I2(int n, double* v, int incv, double scale, int fold, double* sum) {
#elif defined( DSUMI2 )
void dsumI2(int n, double* v, int incv, int fold, double* sum) {
#elif defined( DDOTI2 )
void ddotI2(int n, double* v, int incv, double* y, int incy, int fold, double* sum) {
#endif

	SET_DAZ_FLAG

#ifdef MANUAL_UNROLL
	if (fold == 3 && incv == 1) {
#if defined( DASUMI2 )
		dasumI2_k3(n, v, incv, sum);
		RESET_DAZ_FLAG
		return;
#elif defined( DNRM2I2 )
		dnrm2I2_k3(n, v, incv, scale, sum);
		RESET_DAZ_FLAG
		return;
#elif defined( DSUMI2 )
		dsumI2_k3(n, v, incv, sum);
		RESET_DAZ_FLAG
		return;
#elif defined( DDOTI2 )
		if (incy == 1) {
			ddotI2_k3(n, v, incv, y, incy, sum);
			RESET_DAZ_FLAG
			return;
		}
#endif
	}
#endif
	
	int i, j;
	__m128d mBUF[MAX_FOLD];
	double tmp[2] __attribute__((aligned(16)));
	__m128d mv0, mv1;
	__m128d mq0, mq1;
	__m128d mS0, mR;
	__m128d mBLP;

	// ABSOLUTE MASK
#if defined( DASUMI2 )
	__m128d mAbsMask;
	SSE_ABS_MASKD(mAbsMask);
#elif defined( DNRM2I2 )
	__m128d mScale;
	mScale = _mm_set1_pd(scale);
#elif defined( DDOTI2 )
	__m128d my0, my1;
#endif

	// SET LAST BIT MASK: mBLP = 0.000...1 * 2^0
	SSE_BLP_MASKD(mBLP);

	// EXPAND INITIAL SUM TO BUFFER
	for (j = 0; j < fold; j++) {
		mBUF[j] = _mm_load1_pd(sum + j);
	}

	i = 0;

#if defined( DDOTI2 )
	if (incv == 1 && incy == 1) {
#else
	if (incv == 1) {
#endif
	if (IS_UNALIGNED(v, 16)) {
		mv0 = _mm_load_sd(v);
#if   defined( DASUMI2 )
		mv0 = _mm_and_pd(mv0, mAbsMask);
#elif defined( DNRM2I2 )
		mv0 = _mm_mul_sd(mv0, mScale);
		mv0 = _mm_mul_sd(mv0, mv0);
#elif defined( DDOTI2 )
		my0 = _mm_load_sd(y);
		mv0 = _mm_mul_pd(mv0, my0);
#endif

		for (j = 0; j < fold - 1; j++) {
			mq1     = _mm_or_pd(mv0, mBLP);
			mR     = mBUF[j];
			mBUF[j] = _mm_add_pd(mBUF[j], mq1);
			mR     = _mm_sub_pd(mR, mBUF[j]);
			mv0     = _mm_add_pd(mv0, mR);
		}

		mv0     = _mm_or_pd(mv0, mBLP);
		mBUF[j] = _mm_add_pd(mBUF[j], mv0);
		i = 1;
		v++;
#if defined( DDOTI2 )
		y++;
#endif
	}


#if (DSUMI2_UNROLL_STEP >= 4)
	for (; i < n-3; i+=4, v += 4) {
		mv0 = _mm_load_pd(v);
		mv1 = _mm_load_pd(v + 2);

#if defined( DASUMI2 )
		mv0 = _mm_and_pd(mv0, mAbsMask);
		mv1 = _mm_and_pd(mv1, mAbsMask);
#elif defined( DNRM2I2 )
		mv0 = _mm_mul_pd(mv0, mScale);
		mv0 = _mm_mul_pd(mv0, mv0);
		mv1 = _mm_mul_pd(mv1, mScale);
		mv1 = _mm_mul_pd(mv1, mv1);
#elif defined( DDOTI2 )
		my0 = _mm_loadu_pd(y);
		my1 = _mm_loadu_pd(y + 2);
		mv0 = _mm_mul_pd(mv0, my0);
		mv1 = _mm_mul_pd(mv1, my1);
		y += 4;
#endif

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
	for (; i < n-1; i+=2, v += 2) {
		mv0 = _mm_load_pd(v);
#		if defined( DASUMI2 )
		mv0 = _mm_and_pd(mv0, mAbsMask);
#		elif defined( DNRM2I2 )
		mv0 = _mm_mul_pd(mv0, mScale);
		mv0 = _mm_mul_pd(mv0, mv0);
#		elif defined( DDOTI2 )
		my0 = _mm_loadu_pd(y);
		mv0 = _mm_mul_pd(mv0, my0);
		y += 2;
#endif

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
	}
	else {
	for (; i < n-1; i+=2, v += 2 * incv) {
		mv0 = _mm_load_sd(v);
		mv0 = _mm_loadh_pd(mv0, v + incv);

#		if defined( DASUMI2 )
		mv0 = _mm_and_pd(mv0, mAbsMask);
#		elif defined( DNRM2I2 )
		mv0 = _mm_mul_pd(mv0, mScale);
		mv0 = _mm_mul_pd(mv0, mv0);
#		elif defined( DDOTI2 )
		my0 = _mm_load_sd(y);
		my0 = _mm_loadh_pd(my0, y + incy);
		mv0 = _mm_mul_pd(mv0, my0);
		y += 2 * incy;
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
	}

	// last one
	for (; i < n; i++, v+=incv) {
		mv0 = _mm_load_sd(v);
#		if defined( DASUMI2 )
		mv0 = _mm_and_pd(mv0, mAbsMask);
#		elif defined( DNRM2I2 )
		mv0 = _mm_mul_sd(mv0, mScale);
		mv0 = _mm_mul_sd(mv0, mv0);
#		elif defined( DDOTI2 )
		my0 = _mm_load_sd(y);
		mv0 = _mm_mul_sd(mv0, my0);
		y += incy;
#		endif

		for (j = 0; j < fold - 1; j++) {
			mq1     = _mm_or_pd(mv0, mBLP);
			mq1 = mv0;
			mR     = mBUF[j];
			mBUF[j] = _mm_add_pd(mBUF[j], mq1);
			mR     = _mm_sub_pd(mR, mBUF[j]);
			mv0     = _mm_add_pd(mv0, mR);
		}

		mv0     = _mm_or_pd(mv0, mBLP);
		mBUF[j] = _mm_add_pd(mBUF[j], mv0);
	}


	// TODO SSE3 horizontal add _mm_hadd_pd
	for (j = 0; j < fold; j++) {
		mR = _mm_load_sd(sum + j);
		mBUF[j] = _mm_sub_pd(mBUF[j], mR);
//#ifdef __SSE3__
//		mBUF[j] = _mm_hadd_pd(mBUF[j], mBUF[j]);
//		_mm_store_sd(sum + j, mBUF[j]);
//#else
		_mm_store_pd(tmp, mBUF[j]);
		sum[j] = tmp[0] + tmp[1];
//#endif
	}
	RESET_DAZ_FLAG
}
#endif
