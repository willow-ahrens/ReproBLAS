/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "../types.h"
#include "indexedBLAS.h"

int csum_exception(
	float complex amax,
	int N, float complex* v, int incv,
	int fold, float complex* sum
) {
	int f;
	int sgn = 0;

	// HANDLE NAN AND INIFITY
	sgn = 0;
	if (isinf(CREAL_(amax))) {
		// CHECK FOR NAN
		for (f = 0; f < N; f++) {
			if (isnan(CREAL_(v[f]))) {
				CSET_REAL_(amax, NAN);
				break;
			}
			if (isinf(CREAL_(v[f]))) {
				// POSITIVE INFINITY
				if (signbit(CREAL_(v[f])) == 0) {
					if (sgn < 0) { // INF - INF = NAN
						CSET_REAL_(amax, NAN);
						break;
					}
					sgn = 1;
				}
				else {
					if (sgn > 0) { // INF - INF = NAN
						CSET_REAL_(amax, NAN);
						break;
					}
					sgn = -1;
					CSET_REAL_(amax, -INFINITY);
				}
			}
		}
	}
	sgn = 0;
	if (isinf(CIMAG_(amax))) {
		// CHECK FOR NAN
		for (f = 0; f < N; f++) {
			if (isnan(CIMAG_(v[f]))) {
				CSET_IMAG_(amax, NAN);
				break;
			}
			if (isinf(CIMAG_(v[f]))) {
				// POSITIVE INFINITY
				if (signbit(CIMAG_(v[f])) == 0) {
					if (sgn < 0) { // INF - INF = NAN
						CSET_IMAG_(amax, NAN);
						break;
					}
					sgn = 1;
				}
				else {
					if (sgn > 0) { // INF - INF = NAN
						CSET_IMAG_(amax, NAN);
						break;
					}
					sgn = -1;
					CSET_IMAG_(amax, -INFINITY);
				}
			}
		}
	}

	// NAN
	if (isnan(CREAL_(amax))) {
		for (f = 0; f < fold; f++) 
			CSET_REAL_(sum[f], CREAL_(amax));
	}
	if (isnan(CIMAG_(amax))) {
		for (f = 0; f < fold; f++) 
			CSET_IMAG_(sum[f], CIMAG_(amax));
	}

	// INFINITY
	if (isinf(CREAL_(amax))) {
		if (!isinf(CREAL_(sum[0]))) {
			for (f = 0; f < fold; f++) 
				CSET_REAL_(sum[f], CREAL_(amax));
		}
		if (signbit(CREAL_(amax)) != signbit(CREAL_(sum[0]))) {
			// INF - INF = NAN
			for (f = 0; f < fold; f++) 
				CSET_REAL_(sum[f], NAN);
		}
	}
	// INFINITY
	if (isinf(CIMAG_(amax))) {
		if (!isinf(CIMAG_(sum[0]))) {
			for (f = 0; f < fold; f++) 
				CSET_IMAG_(sum[f], CIMAG_(amax));
		}
		if (signbit(CIMAG_(amax)) != signbit(CIMAG_(sum[0]))) {
			// INF - INF = NAN
			for (f = 0; f < fold; f++) 
				CSET_IMAG_(sum[f], NAN);
		}
	}

	if (isnan(CREAL_(sum[0])) && isnan(CIMAG_(sum[0])))
		return -1;

	if (isinf(CREAL_(sum[0])) && isinf(CIMAG_(sum[0])))
		return 1;

	if (isnan(CREAL_(sum[0])) && isinf(CIMAG_(sum[0])))
		return 1;

	if (isinf(CREAL_(sum[0])) && isnan(CIMAG_(sum[0])))
		return 1;

	if (isnan(CREAL_(sum[0])))
		return 2;

	if (isnan(CIMAG_(sum[0])))
		return 3;

	return 0;
}

void csumI1_(int N, int NB,
	float complex* v, int inc,
	int fold, float complex* sum, float* C) {

	float complex amax;

	int i, j;
	int status;
	int lN;
	int maxN = sICapacity();
	int lB;
		
	for (i = 0; i < N; i+=NB, v+=NB*inc) {
		lN = NB < (N - i) ? NB:(N-i);

		// LOCAL MAX ABSOLUTE
		amax = camax(lN, v, inc);
		status = csum_exception(amax, N, v, inc, fold, sum);

		// TODO: CHECK STATUS
		if (status < 0)
			return;
		if (status == 1)
			continue;
		
		cmcupdate(&amax, sum, 1, C, 1, fold);

		for (j = 0; j < lN; j+=maxN) {
			lB = maxN < lN - j ? maxN : lN - j;
			csumI2(lB, v+j*inc, inc, fold, sum);

			cmrenorm(sum, 1, C, 1, fold);
		}

		if (status == 2) {
			for (j = 0; j < fold; j++)
				CSET_REAL_(sum[j], CREAL_(amax));
		}
		if (status == 3) {
			for (j = 0; j < fold; j++)
				CSET_IMAG_(sum[j], CIMAG_(amax));
		}
	}
}

