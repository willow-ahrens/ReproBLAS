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
#include "IndexedFP/dIndexed.h"
#include "rblas1.h"

#define MANUAL_UNROLL

#define UNROLL_STEP_NR_k3 2

/*
 * Reference SISD implementation for k-fold binning sum/asum/nrm2
 */
#if defined( DZASUMI2 )
void dzasumI2(int n, dcomplex* v, int incv, int fold, dcomplex* sum) {
#endif
#if defined( DZNRM2I2 )
void dznrm2I2(int n, dcomplex* v, int incv, double scale,
	int fold, dcomplex* sum) {
#endif
#if defined( ZSUMI2 )
void zsumI2(int n, dcomplex* v, int incv, int fold, dcomplex* sum) {
#endif

	SET_DAZ_FLAG

	dcomplex zBUF[2 * MAX_FOLD];
	int i, j;

	// EXPAND INITIAL SUM TO BUFFER
	for (j = 0; j < fold; j++) {
		zBUF[2*j] = sum[j];
		zBUF[2*j+1] = sum[j];
	}

	i = 0;
	double* vptr = (double*) v;
	double* sptr = (double*) sum;
	double* BUF  = (double*) zBUF;
	double R0, I0, R1, I1;
	double SR0, SI0;
	double SR1, SI1;
	double SR2, SI2;
	double qR0, qI0;
	l_double ldR0, ldI0;
	l_double ldR1, ldI1;
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
#		if defined( ZSUMI2 )
		R0 = vptr[0];
		I0 = vptr[1];

		R1 = vptr[incv];
		I1 = vptr[incv+1];
#		endif
#		if defined( DZASUMI2 )
		R0 = fabs(vptr[0]);
		I0 = fabs(vptr[1]);

		R1 = fabs(vptr[incv]);
		I1 = fabs(vptr[incv+1]);
#		endif
#		ifdef DZNRM2I2
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
		ldR0.d = R0;
		ldI0.d = I0;
		ldR1.d = R1;
		ldI1.d = I1;

		ldR0.l |= 1;
		ldI0.l |= 1;
		ldR1.l |= 1;
		ldI1.l |= 1;

		qR0 = SR0;
		qI0 = SI0;

		SR0 += ldR0.d;
		SI0 += ldI0.d;

		qR0 -= SR0;
		qI0 -= SI0;

		R0  += qR0;
		I0  += qI0;

		qR0 = SR0;
		qI0 = SI0;

		SR0 += ldR1.d;
		SI0 += ldI1.d;

		qR0 -= SR0;
		qI0 -= SI0;

		R1  += qR0;
		I1  += qI0;
		
		ldR0.d = R0;
		ldI0.d = I0;
		ldR1.d = R1;
		ldI1.d = I1;

		ldR0.l |= 1;
		ldI0.l |= 1;
		ldR1.l |= 1;
		ldI1.l |= 1;

		qR0 = SR1;
		qI0 = SI1;

		SR1 += ldR0.d;
		SI1 += ldI0.d;

		qR0 -= SR1;
		qI0 -= SI1;

		R0  += qR0;
		I0  += qI0;

		qR0 = SR1;
		qI0 = SI1;

		SR1 += ldR1.d;
		SI1 += ldI1.d;

		qR0 -= SR1;
		qI0 -= SI1;

		R1  += qR0;
		I1  += qI0;

		ldR0.d = R0;
		ldI0.d = I0;
		ldR1.d = R1;
		ldI1.d = I1;

		ldR0.l |= 1;
		ldI0.l |= 1;
		ldR1.l |= 1;
		ldI1.l |= 1;

		SR2 += ldR0.d;
		SI2 += ldI0.d;
		SR2 += ldR1.d;
		SI2 += ldI1.d;
	}

	for (; i < n; i++, vptr+=incv) {
#		if defined( ZSUMI2 )
		R0 = vptr[0];
		I0 = vptr[1];
#		endif
#		if defined( DZASUMI2 )
		R0 = fabs(vptr[0]);
		I0 = fabs(vptr[1]);
#		endif
#		ifdef DZNRM2I2
		R0 = scale * vptr[0];
		I0 = scale * vptr[1];

		R0 = R0 * R0;
		I0 = I0 * I0;
#		endif

		// BITWISE OR TO SET THE LAST BIT
		ldR0.d = R0;
		ldI0.d = I0;

		ldR0.l |= 1;
		ldI0.l |= 1;

		qR0 = SR0;
		qI0 = SI0;

		SR0 += ldR0.d;
		SI0 += ldI0.d;

		qR0 -= SR0;
		qI0 -= SI0;

		R0  += qR0;
		I0  += qI0;

		ldR0.d = R0;
		ldI0.d = I0;

		ldR0.l |= 1;
		ldI0.l |= 1;

		qR0 = SR1;
		qI0 = SI1;

		SR1 += ldR0.d;
		SI1 += ldI0.d;

		qR0 -= SR1;
		qI0 -= SI1;

		R0  += qR0;
		I0  += qI0;

		ldR0.d = R0;
		ldI0.d = I0;
		ldR0.l |= 1;
		ldI0.l |= 1;

		SR2 += ldR0.d;
		SI2 += ldI0.d;
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
#		if defined( ZSUMI2 )
		R0 = vptr[0];
		I0 = vptr[1];

		R1 = vptr[incv];
		I1 = vptr[incv+1];
#		endif
#		if defined( DZASUMI2 )
		R0 = fabs(vptr[0]);
		I0 = fabs(vptr[1]);

		R1 = fabs(vptr[incv]);
		I1 = fabs(vptr[incv+1]);
#		endif
#		ifdef DZNRM2I2
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
			ldR0.d = R0;
			ldI0.d = I0;
			ldR0.l |= 1;
			ldI0.l |= 1;

			ldR1.d = R1;
			ldI1.d = I1;
			ldR1.l |= 1;
			ldI1.l |= 1;

			SR0 = sptr[2 * j];
			SI0 = sptr[2 * j + 1];

			qR0 = SR0;
			qI0 = SI0;

			SR0 += ldR0.d;
			SI0 += ldI0.d;

			qR0 -= SR0;
			qI0 -= SI0;

			R0  += qR0;
			I0  += qI0;


			qR0 = SR0;
			qI0 = SI0;

			SR0 += ldR1.d;
			SI0 += ldI1.d;

			qR0 -= SR0;
			qI0 -= SI0;

			R1  += qR0;
			I1  += qI0;

			sptr[2 * j]     = SR0;
			sptr[2 * j + 1] = SI0;
		}

		ldR0.d = R0;
		ldI0.d = I0;
		ldR0.l |= 1;
		ldI0.l |= 1;

		ldR1.d = R1;
		ldI1.d = I1;
		ldR1.l |= 1;
		ldI1.l |= 1;

		sptr[2 * fold - 2] += ldR0.d;
		sptr[2 * fold - 1] += ldI0.d;

		sptr[2 * fold - 2] += ldR1.d;
		sptr[2 * fold - 1] += ldI1.d;

	}

	for (; i < n; i++, vptr+=incv) {
#		if defined( ZSUMI2 )
		R0 = vptr[0];
		I0 = vptr[1];
#		endif
#		if defined( DZASUMI2 )
		R0 = fabs(vptr[0]);
		I0 = fabs(vptr[1]);
#		endif
#		ifdef DZNRM2I2
		R0 = scale * vptr[0];
		I0 = scale * vptr[1];

		R0 = R0 * R0;
		I0 = I0 * I0;
#		endif

		for (j = 0; j < fold - 1; j++) {
			// BITWISE OR TO SET THE LAST BIT
			ldR0.d = R0;
			ldI0.d = I0;
			ldR0.l |= 1;
			ldI0.l |= 1;

			SR0 = sptr[2 * j];
			SI0 = sptr[2 * j + 1];

			qR0 = SR0;
			qI0 = SI0;

			SR0 += ldR0.d;
			SI0 += ldI0.d;

			qR0 -= SR0;
			qI0 -= SI0;

			R0  += qR0;
			I0  += qI0;

			sptr[2 * j]     = SR0;
			sptr[2 * j + 1] = SI0;
		}

		ldR0.d = R0;
		ldI0.d = I0;
		ldR0.l |= 1;
		ldI0.l |= 1;

		sptr[2 * fold - 2] += ldR0.d;
		sptr[2 * fold - 1] += ldI0.d;

	}

	RESET_DAZ_FLAG
}

