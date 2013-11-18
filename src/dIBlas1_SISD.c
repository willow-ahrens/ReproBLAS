/*
 *  Created   13/10/25   H.D. Nguyen
 */

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

#ifndef UNROLL_STEP_RBLAS1
#define UNROLL_STEP_RBLAS1 4
#endif

/*
 * Reference SISD implementation for k-fold binning sum/asum/nrm2/dot
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

	int i, j;
	double BUF[UNROLL_STEP_RBLAS1 * MAX_FOLD];
	double v0, v1;
	double q0, q1;
	double S0, S1, S2;
	double mBLP;
	l_double ld0, ld1;

#ifdef DDOTI2
	double y0, y1;
#endif

#if (UNROLL_STEP_RBLAS1 >= 4)
	double v2, v3;
	double q2, q3;
	double S3;
	l_double ld2, ld3;
	int k;

#ifdef DDOTI2
	double y2, y3;
#endif
#endif

	SET_DAZ_FLAG

	i = 0;

#if defined( DDOTI2 )
	if (incv == 1 && incy == 1) {
#else
	if (incv == 1) {
#endif

#if (UNROLL_STEP_RBLAS1 >= 4)
	if (fold == 3) {
		double S0_1, S0_2, S0_3;
		double S1_1, S1_2, S1_3;
		double S2_1, S2_2, S2_3;

		S0_1 = S0_2 = S0_3 = S0 = sum[0];
		S1_1 = S1_2 = S1_3 = S1 = sum[1];
		S2_1 = S2_2 = S2_3 = S2 = sum[2];
		for (; i < n-3; i+=4, v += 4) {

#if defined( DSUMI2 )
			v0 = (v[0]);
			v1 = (v[1]);
			v2 = (v[2]);
			v3 = (v[3]);
#elif defined( DASUMI2 )
			v0 = fabs(v[0]);
			v1 = fabs(v[1]);
			v2 = fabs(v[2]);
			v3 = fabs(v[3]);
#elif defined( DNRM2I2 )
			v0 = v[0] * scale;
			v0 *= v0;

			v1 =  v[1] * scale;
			v1 *= v1;

			v2 =  v[2] * scale;
			v2 *= v2;

			v3  = v[3] * scale;
			v3 *= v3;
#elif defined( DDOTI2 )
			v0 = v[0] * y[0];
			v1 = v[1] * y[1];
			v2 = v[2] * y[2];
			v3 = v[3] * y[3];

			y += 4;
#endif

			// set last bit
			ld0.d  = v0;
			ld0.l |= 1;
			//...
			ld1.d  = v1;
			ld1.l |= 1;
			//...
			ld2.d  = v2;
			ld2.l |= 1;
			//...
			ld3.d  = v3;
			ld3.l |= 1;
			//----

			q0     = S0;
			S0     = S0 + ld0.d;
			q0     = q0 - S0;
			v0     = v0 + q0;
			//...
			q1     = S0_1;
			S0_1   = S0_1 + ld1.d;
			q1     = q1 - S0_1;
			v1     = v1 + q1;
			//...
			q2     = S0_2;
			S0_2    = S0_2 + ld2.d;
			q2     = q2 - S0_2;
			v2     = v2 + q2;
			//...
			q3     = S0_3;
			S0_3   = S0_3 + ld3.d;
			q3     = q3 - S0_3;
			v3     = v3 + q3;
			//---

			// set last bit
			ld0.d  = v0;
			ld0.l |= 1;
			//...
			ld1.d  = v1;
			ld1.l |= 1;
			//...
			ld2.d  = v2;
			ld2.l |= 1;
			//...
			ld3.d  = v3;
			ld3.l |= 1;
			//----

			q0     = S1;
			S1     = S1 + ld0.d;
			q0     = q0 - S1;
			v0     = v0 + q0;
			//...
			q1     = S1_1;
			S1_1   = S1_1 + ld1.d;
			q1     = q1 - S1_1;
			v1     = v1 + q1;
			//...
			q2     = S1_2;
			S1_2    = S1_2 + ld2.d;
			q2     = q2 - S1_2;
			v2     = v2 + q2;
			//...
			q3     = S1_3;
			S1_3   = S1_3 + ld3.d;
			q3     = q3 - S1_3;
			v3     = v3 + q3;
			//---

			// set last bit
			ld0.d  = v0;
			ld0.l |= 1;
			//...
			ld1.d  = v1;
			ld1.l |= 1;
			//...
			ld2.d  = v2;
			ld2.l |= 1;
			//...
			ld3.d  = v3;
			ld3.l |= 1;
			//----

			S2   += ld0.d;
			S2_1 += ld1.d;
			S2_2 += ld2.d;
			S2_3 += ld3.d;
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

#if defined( DSUMI2 )
		v0 = (v[0]);
		v1 = (v[1]);
		v2 = (v[2]);
		v3 = (v[3]);
#elif defined( DASUMI2 )
		v0 = fabs(v[0]);
		v1 = fabs(v[1]);
		v2 = fabs(v[2]);
		v3 = fabs(v[3]);
#elif defined( DNRM2I2 )
		v0 = v[0] * scale;
		v0 *= v0;

		v1 =  v[1] * scale;
		v1 *= v1;

		v2 =  v[2] * scale;
		v2 *= v2;

		v3  = v[3] * scale;
		v3 *= v3;
#elif defined( DDOTI2 )
		v0 = v[0] * y[0];
		v1 = v[1] * y[1];
		v2 = v[2] * y[2];
		v3 = v[3] * y[3];

		y += 4;
#endif

		for (j = 0; j < fold - 1; j++) {
			// set last bit
			ld0.d  = v0;
			ld0.l |= 1;
			q0     = ld0.d;
			//...
			ld1.d  = v1;
			ld1.l |= 1;
			q1     = ld1.d;
			//...
			ld2.d  = v2;
			ld2.l |= 1;
			q2     = ld2.d;
			//...
			ld3.d  = v3;
			ld3.l |= 1;
			q3     = ld3.d;
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
		ld0.d  = v0;
		ld0.l |= 1;
		q0     = ld0.d;
		//...
		ld1.d  = v1;
		ld1.l |= 1;
		q1     = ld1.d;
		//...
		ld2.d  = v2;
		ld2.l |= 1;
		q2     = ld2.d;
		//...
		ld3.d  = v3;
		ld3.l |= 1;
		q3     = ld3.d;
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
#       	if   defined (DSUMI2)
			v0 = (v[0]);
#			elif defined( DASUMI2 )
			v0 = fabs(v[0]);
#			elif defined( DNRM2I2 )
			v0 = v[0] * scale;
			v0 = v0 * v0;
#			elif defined( DDOTI2 )
			v0 = v[0] * y[0];
			y += incy;
#			endif

			// set last bit
			ld1.d  = v0;
			ld1.l |= 1;

			S0 += ld1.d;
		}
		sum[0] = S0;
	}
	// 2 BINS
	else if (fold == 2) {
		S0 = sum[0];
		S1 = sum[1];
		for (; i < n; i++, v+=incv) {
#       	if   defined (DSUMI2)
			v0 = (v[0]);
#			elif defined( DASUMI2 )
			v0 = fabs(v[0]);
#			elif defined( DNRM2I2 )
			v0 = v[0] * scale;
			v0 = v0 * v0;
#			elif defined( DDOTI2 )
			v0 = v[0] * y[0];
			y += incy;
#			endif

			// set last bit
			ld1.d  = v0;
			ld1.l |= 1;

			q0     = S0;
			S0     = S0 + ld1.d;
			q0     = q0 - S0;
			v0     = v0 + S0;

			// set last bit
			ld1.d  = v0;
			ld1.l |= 1;

			S1 += ld1.d;
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
#       	if   defined (DSUMI2)
			v0 = (v[0]);
#			elif defined( DASUMI2 )
			v0 = fabs(v[0]);
#			elif defined( DNRM2I2 )
			v0 = v[0] * scale;
			v0 = v0 * v0;
#			elif defined( DDOTI2 )
			v0 = v[0] * y[0];
			y += incy;
#			endif

			// set last bit
			ld1.d  = v0;
			ld1.l |= 1;

			q0     = S0;
			S0     = S0 + ld1.d;
			q0     = q0 - S0;
			v0     = v0 + q0;

			// set last bit
			ld1.d  = v0;
			ld1.l |= 1;

			q0     = S1;
			S1     = S1 + ld1.d;
			q0     = q0 - S1;
			v0     = v0 + q0;

			// set last bit
			ld1.d  = v0;
			ld1.l |= 1;

			S2 += ld1.d;
		}
		sum[0] = S0;
		sum[1] = S1;
		sum[2] = S2;
	}
	else
	for (; i < n; i++, v+=incv) {
#       if   defined (DSUMI2)
		v0 = (v[0]);
#		elif defined( DASUMI2 )
		v0 = fabs(v[0]);
#		elif defined( DNRM2I2 )
		v0 = v[0] * scale;
		v0 = v0 * v0;
#		elif defined( DDOTI2 )
		v0 = v[0] * y[0];
		y += incy;
#		endif

		for (j = 0; j < fold - 1; j++) {
			// set last bit
			ld1.d  = v0;
			ld1.l |= 1;
			q0     = ld1.d;

			S0     = sum[j];
			q0     = S0 + q0;
			sum[j] = q0;
			S0     = S0 - q0;
			v0     = v0 + S0;
		}
		// set last bit
		ld1.d  = v0;
		ld1.l |= 1;
		q0     = ld1.d;

		sum[j] += q0;
	}

	RESET_DAZ_FLAG

}

