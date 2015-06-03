/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <emmintrin.h>

#include "../config.h"
#include "../common/common.h"
#include "cIBlas1.h"
//#include "../indexedFP/sindexed.h"
#include "../rblas1.h"

#define MANUAL_UNROLL

#define UNROLL_STEP_NR_k3 2

/*
 * Reference SISD implementation for k-fold binning sum/asum/nrm2
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

	float complex zBUF[2 * MAX_FOLD];
	int i, j;

	// EXPAND INITIAL SUM TO BUFFER
	for (j = 0; j < fold; j++) {
		zBUF[2*j] = sum[j];
		zBUF[2*j+1] = sum[j];
	}

	i = 0;
	float* vptr = (float*) v;
	float* sptr = (float*) sum;
	float* BUF  = (float*) zBUF;
	float R0, I0, R1, I1;
	float SR0, SI0;
	float SR1, SI1;
	float SR2, SI2;
	float qR0, qI0;
	int_float ldR0, ldI0;
	int_float ldR1, ldI1;
	incv *= 2;

	//++++ 3-fold ++++
	if (fold == 3) {
		SR0 = sptr[0];
		SI0 = sptr[1];
		SR1 = sptr[2];
		SI1 = sptr[3];
		SR2 = sptr[4];
		SI2 = sptr[5];

	for (; i < n-1; i+=2, vptr+=2*incv) {
#		if defined( CSUMI2 )
		R0 = vptr[0];
		I0 = vptr[1];

		R1 = vptr[incv];
		I1 = vptr[incv+1];
#		endif
#		if defined( SCASUMI2 )
		R0 = fabs(vptr[0]);
		I0 = fabs(vptr[1]);

		R1 = fabs(vptr[incv]);
		I1 = fabs(vptr[incv+1]);
#		endif
#		ifdef SCNRM2I2
		R0 = scale * vptr[0];
		I0 = scale * vptr[1];

		R0 = R0 * R0;
		I0 = I0 * I0;

		R1 = scale * vptr[incv];
		I1 = scale * vptr[incv+1];

		R1 = R1 * R1;
		I1 = I1 * I1;
#		endif

		// BITWISE OR TO SET THE LAST BIT
		ldR0.f = R0;
		ldI0.f = I0;
		ldR1.f = R1;
		ldI1.f = I1;

		ldR0.i |= 1;
		ldI0.i |= 1;
		ldR1.i |= 1;
		ldI1.i |= 1;

		qR0 = SR0;
		qI0 = SI0;

		SR0 += ldR0.f;
		SI0 += ldI0.f;

		qR0 -= SR0;
		qI0 -= SI0;

		R0  += qR0;
		I0  += qI0;

		qR0 = SR0;
		qI0 = SI0;

		SR0 += ldR1.f;
		SI0 += ldI1.f;

		qR0 -= SR0;
		qI0 -= SI0;

		R1  += qR0;
		I1  += qI0;
		
		ldR0.f = R0;
		ldI0.f = I0;
		ldR1.f = R1;
		ldI1.f = I1;

		ldR0.i |= 1;
		ldI0.i |= 1;
		ldR1.i |= 1;
		ldI1.i |= 1;

		qR0 = SR1;
		qI0 = SI1;

		SR1 += ldR0.f;
		SI1 += ldI0.f;

		qR0 -= SR1;
		qI0 -= SI1;

		R0  += qR0;
		I0  += qI0;

		qR0 = SR1;
		qI0 = SI1;

		SR1 += ldR1.f;
		SI1 += ldI1.f;

		qR0 -= SR1;
		qI0 -= SI1;

		R1  += qR0;
		I1  += qI0;

		ldR0.f = R0;
		ldI0.f = I0;
		ldR1.f = R1;
		ldI1.f = I1;

		ldR0.i |= 1;
		ldI0.i |= 1;
		ldR1.i |= 1;
		ldI1.i |= 1;

		SR2 += ldR0.f;
		SI2 += ldI0.f;
		SR2 += ldR1.f;
		SI2 += ldI1.f;
	}

	for (; i < n; i++, vptr+=incv) {
#		if defined( CSUMI2 )
		R0 = vptr[0];
		I0 = vptr[1];
#		endif
#		if defined( SCASUMI2 )
		R0 = fabs(vptr[0]);
		I0 = fabs(vptr[1]);
#		endif
#		ifdef SCNRM2I2
		R0 = scale * vptr[0];
		I0 = scale * vptr[1];

		R0 = R0 * R0;
		I0 = I0 * I0;
#		endif

		// BITWISE OR TO SET THE LAST BIT
		ldR0.f = R0;
		ldI0.f = I0;

		ldR0.i |= 1;
		ldI0.i |= 1;

		qR0 = SR0;
		qI0 = SI0;

		SR0 += ldR0.f;
		SI0 += ldI0.f;

		qR0 -= SR0;
		qI0 -= SI0;

		R0  += qR0;
		I0  += qI0;

		ldR0.f = R0;
		ldI0.f = I0;

		ldR0.i |= 1;
		ldI0.i |= 1;

		qR0 = SR1;
		qI0 = SI1;

		SR1 += ldR0.f;
		SI1 += ldI0.f;

		qR0 -= SR1;
		qI0 -= SI1;

		R0  += qR0;
		I0  += qI0;

		ldR0.f = R0;
		ldI0.f = I0;
		ldR0.i |= 1;
		ldI0.i |= 1;

		SR2 += ldR0.f;
		SI2 += ldI0.f;
	}

	RESET_DAZ_FLAG

	sptr[0] = SR0;
	sptr[1] = SI0;
	sptr[2] = SR1;
	sptr[3] = SI1;
	sptr[4] = SR2;
	sptr[5] = SI2;
	return;
	}
	//---- 3-fold ----

	for (; i < n-1; i+=2, vptr+=2*incv) {
#		if defined( CSUMI2 )
		R0 = vptr[0];
		I0 = vptr[1];

		R1 = vptr[incv];
		I1 = vptr[incv+1];
#		endif
#		if defined( SCASUMI2 )
		R0 = fabs(vptr[0]);
		I0 = fabs(vptr[1]);

		R1 = fabs(vptr[incv]);
		I1 = fabs(vptr[incv+1]);
#		endif
#		ifdef SCNRM2I2
		R0 = scale * vptr[0];
		I0 = scale * vptr[1];

		R0 = R0 * R0;
		I0 = I0 * I0;

		R1 = scale * vptr[incv];
		I1 = scale * vptr[incv+1];

		R1 = R1 * R1;
		I1 = I1 * I1;
#		endif

		for (j = 0; j < fold - 1; j++) {
			// BITWISE OR TO SET THE LAST BIT
			ldR0.f = R0;
			ldI0.f = I0;
			ldR0.i |= 1;
			ldI0.i |= 1;

			ldR1.f = R1;
			ldI1.f = I1;
			ldR1.i |= 1;
			ldI1.i |= 1;

			SR0 = sptr[2 * j];
			SI0 = sptr[2 * j + 1];

			qR0 = SR0;
			qI0 = SI0;

			SR0 += ldR0.f;
			SI0 += ldI0.f;

			qR0 -= SR0;
			qI0 -= SI0;

			R0  += qR0;
			I0  += qI0;


			qR0 = SR0;
			qI0 = SI0;

			SR0 += ldR1.f;
			SI0 += ldI1.f;

			qR0 -= SR0;
			qI0 -= SI0;

			R1  += qR0;
			I1  += qI0;

			sptr[2 * j]     = SR0;
			sptr[2 * j + 1] = SI0;
		}

		ldR0.f = R0;
		ldI0.f = I0;
		ldR0.i |= 1;
		ldI0.i |= 1;

		ldR1.f = R1;
		ldI1.f = I1;
		ldR1.i |= 1;
		ldI1.i |= 1;

		sptr[2 * fold - 2] += ldR0.f;
		sptr[2 * fold - 1] += ldI0.f;

		sptr[2 * fold - 2] += ldR1.f;
		sptr[2 * fold - 1] += ldI1.f;

	}

	for (; i < n; i++, vptr+=incv) {
#		if defined( CSUMI2 )
		R0 = vptr[0];
		I0 = vptr[1];
#		endif
#		if defined( SCASUMI2 )
		R0 = fabs(vptr[0]);
		I0 = fabs(vptr[1]);
#		endif
#		ifdef SCNRM2I2
		R0 = scale * vptr[0];
		I0 = scale * vptr[1];

		R0 = R0 * R0;
		I0 = I0 * I0;
#		endif

		for (j = 0; j < fold - 1; j++) {
			// BITWISE OR TO SET THE LAST BIT
			ldR0.f = R0;
			ldI0.f = I0;
			ldR0.i |= 1;
			ldI0.i |= 1;

			SR0 = sptr[2 * j];
			SI0 = sptr[2 * j + 1];

			qR0 = SR0;
			qI0 = SI0;

			SR0 += ldR0.f;
			SI0 += ldI0.f;

			qR0 -= SR0;
			qI0 -= SI0;

			R0  += qR0;
			I0  += qI0;

			sptr[2 * j]     = SR0;
			sptr[2 * j + 1] = SI0;
		}

		ldR0.f = R0;
		ldI0.f = I0;
		ldR0.i |= 1;
		ldI0.i |= 1;

		sptr[2 * fold - 2] += ldR0.f;
		sptr[2 * fold - 1] += ldI0.f;

	}

	RESET_DAZ_FLAG
}

