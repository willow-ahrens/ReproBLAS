/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <emmintrin.h>

#include "../config.h"
#include "../Common/Common.h"
#include "zIBlas1.h"
#include "../types.h"

#define MANUAL_UNROLL

#define UNROLL_STEP_NR_k3 2

#define DSET_LAST_BIT(X,XB) \
{	\
	long_double ld;	\
	ld.d  = X;		\
	ld.l |= 1;		\
	XB    = ld.d;	\
}

/*
 * Reference SISD implementation for k-fold binning sum/asum/nrm2
 */
#ifdef ZDOTCI2
void zdotcI2
#elif defined (ZDOTUI2)
void zdotuI2
#endif
(int n, dcomplex* v, int incv, dcomplex* y, int incy, int fold, dcomplex* sum)
{
	SET_DAZ_FLAG

	double RR[2 * MAX_FOLD];
	double II[2 * MAX_FOLD];
	double R0, I0, R1, I1;
	double RR0, II0, R0I1, R1I0;
	double SR0, SR1, SI0, SI1;
	double real, imag;

	int i, j;
	double* vptr = (double*)v;
	double* yptr = (double*)y;
	dcomplex* zRR = (dcomplex*) RR;
	dcomplex* zII = (dcomplex*) II;

	// EXPAND INITIAL SUM TO BUFFER
	for (j = 0; j < fold; j++) {
		RR[2*j]   = ZREAL_(sum[j]);
		RR[2*j+1] = RR[2*j];
		II[2*j]   = ZIMAG_(sum[j]);
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
			DSET_LAST_BIT(RR0, R0);
			DSET_LAST_BIT(II0, I0);
			DSET_LAST_BIT(R1I0, R1);
			DSET_LAST_BIT(R0I1, I1);

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

		DSET_LAST_BIT(RR0, R0);
		DSET_LAST_BIT(II0, I0);
		DSET_LAST_BIT(R1I0, R1);
		DSET_LAST_BIT(R0I1, I1);

		RR[2*j]   += R0;
		RR[2*j+1] += I0;
		II[2*j]   += R1;
		II[2*j+1] += I1;
	}
	
	for (j = 0; j < fold; j++) {

#ifdef ZDOTUI2
		real = (ZREAL_(sum[j]) - ZIMAG_(zRR[j])) + ZREAL_(zRR[j]);
#elif defined ( ZDOTCI2 )
		real = (ZIMAG_(zRR[j]) - ZREAL_(sum[j])) + ZREAL_(zRR[j]);
#endif

#ifdef ZDOTCI2
		imag = (ZIMAG_(sum[j]) - ZIMAG_(zII[j])) + ZREAL_(zII[j]);
#elif defined ( ZDOTUI2 )
		imag = (ZIMAG_(zII[j]) - ZIMAG_(sum[j])) + ZREAL_(zII[j]);
#endif

		ZSET_(sum[j], real, imag);
	}

	RESET_DAZ_FLAG
}

