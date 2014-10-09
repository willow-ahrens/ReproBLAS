/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <emmintrin.h>

#include "../config.h"
#include "../Common/Common.h"
#include "sIBlas1.h"

//#define SASUMI2

/*
 * SUMMATION IN INDEXEDFP FORMAT
 * REFERENCE IMPLEMENTATION
 */

#ifndef UNROLL_STEP_RSBLAS1
#define UNROLL_STEP_RSBLAS1 4
#endif

/*
 * Reference SISD implementation for k-fold binning sum/asum/nrm2/dot
 */
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
	float BUF[UNROLL_STEP_RSBLAS1 * MAX_FOLD];
	float v0, v1;
	float q0, q1;
	float S0, S1, S2;
	i_float ld0, ld1;

#if (UNROLL_STEP_RSBLAS1 >= 4)
	float v2, v3;
	float q2, q3;
	float S3;
	i_float ld2, ld3;

#endif

	SET_DAZ_FLAG

	i = 0;

#if defined( SDOTI2 )
	if (incv == 1 && incy == 1) {
#else
	if (incv == 1) {
#endif

#if (UNROLL_STEP_RSBLAS1 >= 4)
	if (fold == 3) {
		float S0_1, S0_2, S0_3;
		float S1_1, S1_2, S1_3;
		float S2_1, S2_2, S2_3;

		S0_1 = S0_2 = S0_3 = S0 = sum[0];
		S1_1 = S1_2 = S1_3 = S1 = sum[1];
		S2_1 = S2_2 = S2_3 = S2 = sum[2];
		for (; i < n-3; i+=4, v += 4) {

#if defined( SSUMI2 )
			v0 = (v[0]);
			v1 = (v[1]);
			v2 = (v[2]);
			v3 = (v[3]);
#elif defined( SASUMI2 )
			v0 = fabs(v[0]);
			v1 = fabs(v[1]);
			v2 = fabs(v[2]);
			v3 = fabs(v[3]);
#elif defined( SNRM2I2 )
			v0 = v[0] * scale;
			v0 *= v0;

			v1 =  v[1] * scale;
			v1 *= v1;

			v2 =  v[2] * scale;
			v2 *= v2;

			v3  = v[3] * scale;
			v3 *= v3;
#elif defined( SDOTI2 )
			v0 = v[0] * y[0];
			v1 = v[1] * y[1];
			v2 = v[2] * y[2];
			v3 = v[3] * y[3];

			y += 4;
#endif

			// set last bit
			ld0.f  = v0;
			ld0.i |= 1;
			//...
			ld1.f  = v1;
			ld1.i |= 1;
			//...
			ld2.f  = v2;
			ld2.i |= 1;
			//...
			ld3.f  = v3;
			ld3.i |= 1;
			//----

			q0     = S0;
			S0     = S0 + ld0.f;
			q0     = q0 - S0;
			v0     = v0 + q0;
			//...
			q1     = S0_1;
			S0_1   = S0_1 + ld1.f;
			q1     = q1 - S0_1;
			v1     = v1 + q1;
			//...
			q2     = S0_2;
			S0_2    = S0_2 + ld2.f;
			q2     = q2 - S0_2;
			v2     = v2 + q2;
			//...
			q3     = S0_3;
			S0_3   = S0_3 + ld3.f;
			q3     = q3 - S0_3;
			v3     = v3 + q3;
			//---

			// set last bit
			ld0.f  = v0;
			ld0.i |= 1;
			//...
			ld1.f  = v1;
			ld1.i |= 1;
			//...
			ld2.f  = v2;
			ld2.i |= 1;
			//...
			ld3.f  = v3;
			ld3.i |= 1;
			//----

			q0     = S1;
			S1     = S1 + ld0.f;
			q0     = q0 - S1;
			v0     = v0 + q0;
			//...
			q1     = S1_1;
			S1_1   = S1_1 + ld1.f;
			q1     = q1 - S1_1;
			v1     = v1 + q1;
			//...
			q2     = S1_2;
			S1_2    = S1_2 + ld2.f;
			q2     = q2 - S1_2;
			v2     = v2 + q2;
			//...
			q3     = S1_3;
			S1_3   = S1_3 + ld3.f;
			q3     = q3 - S1_3;
			v3     = v3 + q3;
			//---

			// set last bit
			ld0.f  = v0;
			ld0.i |= 1;
			//...
			ld1.f  = v1;
			ld1.i |= 1;
			//...
			ld2.f  = v2;
			ld2.i |= 1;
			//...
			ld3.f  = v3;
			ld3.i |= 1;
			//----

			S2   += ld0.f;
			S2_1 += ld1.f;
			S2_2 += ld2.f;
			S2_3 += ld3.f;
		}
		S0 += (S0_1 - sum[0]);
		S0 += (S0_2 - sum[0]);
		S0 += (S0_3 - sum[0]);
		sum[0] = S0;
		S1 += (S1_1 - sum[1]);
		S1 += (S1_2 - sum[1]);
		S1 += (S1_3 - sum[1]);
		sum[1] = S1;
		S2 += (S2_1 - sum[2]);
		S2 += (S2_2 - sum[2]);
		S2 += (S2_3 - sum[2]);
		sum[2] = S2;
	}
	else {
	// EXPAND INITIAL SUM TO BUFFER
	for (j = 0; j < fold; j++) {
		BUF[j * 4] = sum[ j ];
		BUF[j * 4+1] = sum[ j ];
		BUF[j * 4+2] = sum[ j ];
		BUF[j * 4+3] = sum[ j ];
	}
	for (; i < n-3; i+=4, v += 4) {

#if defined( SSUMI2 )
		v0 = (v[0]);
		v1 = (v[1]);
		v2 = (v[2]);
		v3 = (v[3]);
#elif defined( SASUMI2 )
		v0 = fabs(v[0]);
		v1 = fabs(v[1]);
		v2 = fabs(v[2]);
		v3 = fabs(v[3]);
#elif defined( SNRM2I2 )
		v0 = v[0] * scale;
		v0 *= v0;

		v1 =  v[1] * scale;
		v1 *= v1;

		v2 =  v[2] * scale;
		v2 *= v2;

		v3  = v[3] * scale;
		v3 *= v3;
#elif defined( SDOTI2 )
		v0 = v[0] * y[0];
		v1 = v[1] * y[1];
		v2 = v[2] * y[2];
		v3 = v[3] * y[3];

		y += 4;
#endif

		for (j = 0; j < fold - 1; j++) {
			// set last bit
			ld0.f  = v0;
			ld0.i |= 1;
			q0     = ld0.f;
			//...
			ld1.f  = v1;
			ld1.i |= 1;
			q1     = ld1.f;
			//...
			ld2.f  = v2;
			ld2.i |= 1;
			q2     = ld2.f;
			//...
			ld3.f  = v3;
			ld3.i |= 1;
			q3     = ld3.f;
			//----

			S0     = BUF[4*j];
			S1     = BUF[4*j + 1];
			S2     = BUF[4*j + 2];
			S3     = BUF[4*j + 3];

			q0     = S0 + q0;
			q1     = S1 + q1;
			q2     = S2 + q2;
			q3     = S3 + q3;

			BUF[4*j]   = q0;
			BUF[4*j+1] = q1;
			BUF[4*j+2] = q2;
			BUF[4*j+3] = q3;

			S0     = S0 - q0;
			S1     = S1 - q1;
			S2     = S2 - q2;
			S3     = S3 - q3;

			v0     = v0 + S0;
			v1     = v1 + S1;
			v2     = v2 + S2;
			v3     = v3 + S3;
		}
		// set last bit
		ld0.f  = v0;
		ld0.i |= 1;
		q0     = ld0.f;
		//...
		ld1.f  = v1;
		ld1.i |= 1;
		q1     = ld1.f;
		//...
		ld2.f  = v2;
		ld2.i |= 1;
		q2     = ld2.f;
		//...
		ld3.f  = v3;
		ld3.i |= 1;
		q3     = ld3.f;
		//----

		BUF[4*j]   += q0;
		BUF[4*j+1] += q1;
		BUF[4*j+2] += q2;
		BUF[4*j+3] += q3;
	}

	for (j = 0; j < fold; j++) {
		S0 = BUF[4*j];
		S0 += (BUF[4*j+1] - sum[j]);
		S0 += (BUF[4*j+2] - sum[j]);
		S0 += (BUF[4*j+3] - sum[j]);
		sum[j] = S0;
	}
	}
#endif

	}

	// last one
	// 1 BIN
	if (fold == 1) {
		S0 = sum[0];
		for (; i < n; i++, v+=incv) {
#       	if   defined (SSUMI2)
			v0 = (v[0]);
#			elif defined( SASUMI2 )
			v0 = fabs(v[0]);
#			elif defined( SNRM2I2 )
			v0 = v[0] * scale;
			v0 = v0 * v0;
#			elif defined( SDOTI2 )
			v0 = v[0] * y[0];
			y += incy;
#			endif

			// set last bit
			ld1.f  = v0;
			ld1.i |= 1;

			S0 += ld1.f;
		}
		sum[0] = S0;
	}
	// 2 BINS
	else if (fold == 2) {
		S0 = sum[0];
		S1 = sum[1];
		for (; i < n; i++, v+=incv) {
#       	if   defined (SSUMI2)
			v0 = (v[0]);
#			elif defined( SASUMI2 )
			v0 = fabs(v[0]);
#			elif defined( SNRM2I2 )
			v0 = v[0] * scale;
			v0 = v0 * v0;
#			elif defined( SDOTI2 )
			v0 = v[0] * y[0];
			y += incy;
#			endif

			// set last bit
			ld1.f  = v0;
			ld1.i |= 1;

			q0     = S0;
			S0     = S0 + ld1.f;
			q0     = q0 - S0;
			v0     = v0 + S0;

			// set last bit
			ld1.f  = v0;
			ld1.i |= 1;

			S1 += ld1.f;
		}
		sum[0] = S0;
		sum[1] = S1;
	}
	// 3 BINS
	else if (fold == 3) {
		S0 = sum[0];
		S1 = sum[1];
		S2 = sum[2];
		for (; i < n; i++, v+=incv) {
#       	if   defined (SSUMI2)
			v0 = (v[0]);
#			elif defined( SASUMI2 )
			v0 = fabs(v[0]);
#			elif defined( SNRM2I2 )
			v0 = v[0] * scale;
			v0 = v0 * v0;
#			elif defined( SDOTI2 )
			v0 = v[0] * y[0];
			y += incy;
#			endif

			// set last bit
			ld1.f  = v0;
			ld1.i |= 1;

			q0     = S0;
			S0     = S0 + ld1.f;
			q0     = q0 - S0;
			v0     = v0 + q0;

			// set last bit
			ld1.f  = v0;
			ld1.i |= 1;

			q0     = S1;
			S1     = S1 + ld1.f;
			q0     = q0 - S1;
			v0     = v0 + q0;

			// set last bit
			ld1.f  = v0;
			ld1.i |= 1;

			S2 += ld1.f;
		}
		sum[0] = S0;
		sum[1] = S1;
		sum[2] = S2;
	}
	else
	for (; i < n; i++, v+=incv) {
#       if   defined (SSUMI2)
		v0 = (v[0]);
#		elif defined( SASUMI2 )
		v0 = fabs(v[0]);
#		elif defined( SNRM2I2 )
		v0 = v[0] * scale;
		v0 = v0 * v0;
#		elif defined( SDOTI2 )
		v0 = v[0] * y[0];
		y += incy;
#		endif

		for (j = 0; j < fold - 1; j++) {
			// set last bit
			ld1.f  = v0;
			ld1.i |= 1;
			q0     = ld1.f;

			S0     = sum[j];
			q0     = S0 + q0;
			sum[j] = q0;
			S0     = S0 - q0;
			v0     = v0 + S0;
		}
		// set last bit
		ld1.f  = v0;
		ld1.i |= 1;
		q0     = ld1.f;

		sum[j] += q0;
	}

	RESET_DAZ_FLAG

}

