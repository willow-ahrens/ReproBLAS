/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <emmintrin.h>

#include "../config.h"
#include "../Common/Common.h"
#include "cIBlas1.h"
#include "../types.h"

#define MANUAL_UNROLL

#define UNROLL_STEP_NR_k3 2

#define CSET_LAST_BIT(X,XB) \
{	\
	int_float ld;	\
	ld.f  = X;		\
	ld.i |= 1;		\
	XB    = ld.f;	\
}

/*
 * Reference SISD implementation for k-fold binning sum/asum/nrm2
 */
#ifdef CDOTCI2
void cdotcI2
#elif defined (CDOTUI2)
void cdotuI2
#endif
(int n, float complex* v, int incv, float complex* y, int incy, int fold, float complex* sum)
{
	SET_DAZ_FLAG

	float RR[2 * MAX_FOLD];
	float II[2 * MAX_FOLD];
	float R0, I0, R1, I1;
	float RR0, II0, R0I1, R1I0;
	float SR0, SR1, SI0, SI1;
	float real, imag;

	int i, j;
	float* vptr = (float*)v;
	float* yptr = (float*)y;
	float complex* zRR = (float complex*) RR;
	float complex* zII = (float complex*) II;

	// EXPAND INITIAL SUM TO BUFFER
	for (j = 0; j < fold; j++) {
		RR[2*j]   = CREAL_(sum[j]);
		RR[2*j+1] = RR[2*j];
		II[2*j]   = CIMAG_(sum[j]);
		II[2*j+1] = II[2*j];
	}

	i = 0;
	incv *= 2;
	incy *= 2;

	for (; i < n; i++, vptr += incv, yptr += incy) {
		R0 = vptr[0];
		I0 = vptr[1];
		R1 = yptr[0];
		I1 = yptr[1];

		RR0 = R0 * R1;
		II0 = I0 * I1;
		R0I1 = R0 * I1;
		R1I0 = R1 * I0;

		for (j = 0; j < fold - 1; j++) {
			CSET_LAST_BIT(RR0, R0);
			CSET_LAST_BIT(II0, I0);
			CSET_LAST_BIT(R1I0, R1);
			CSET_LAST_BIT(R0I1, I1);

			SR0   = RR[2*j];
			SR1   = RR[2*j + 1];
			SI0   = II[2*j];
			SI1   = II[2*j + 1];

			R0   += SR0;
			I0   += SR1;
			R1   += SI0;
			I1   += SI1;

			SR0  -= R0;
			SR1  -= I0;
			SI0  -= R1;
			SI1  -= I1;

			RR0  += SR0;
			II0  += SR1;
			R1I0 += SI0;
			R0I1 += SI1;

			RR[2*j]   = R0;
			RR[2*j+1] = I0;
			II[2*j]   = R1;
			II[2*j+1] = I1;
		}

		CSET_LAST_BIT(RR0, R0);
		CSET_LAST_BIT(II0, I0);
		CSET_LAST_BIT(R1I0, R1);
		CSET_LAST_BIT(R0I1, I1);

		RR[2*j]   += R0;
		RR[2*j+1] += I0;
		II[2*j]   += R1;
		II[2*j+1] += I1;
	}
	
	for (j = 0; j < fold; j++) {

#ifdef CDOTUI2
		real = (CREAL_(sum[j]) - CIMAG_(zRR[j])) + CREAL_(zRR[j]);
#elif defined ( CDOTCI2 )
		real = (CIMAG_(zRR[j]) - CREAL_(sum[j])) + CREAL_(zRR[j]);
#endif

#ifdef CDOTCI2
		imag = (CIMAG_(sum[j]) - CIMAG_(zII[j])) + CREAL_(zII[j]);
#elif defined ( CDOTUI2 )
		imag = (CIMAG_(zII[j]) - CIMAG_(sum[j])) + CREAL_(zII[j]);
#endif

		CSET_(sum[j], real, imag);
	}

	RESET_DAZ_FLAG
}

