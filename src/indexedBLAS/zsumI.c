/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexedBLAS.h"
#include "../types.h"

// TODO CHECK FOR FALSE INFINITY WHEN ABSOLUTE MAX IS CLOSE TO INFINITY
// SOLUTION: SCALE INPUTS

int zsum_exception(
	double complex amax,
	int N, double complex* v, int incv,
	int fold, double complex* sum
) {
	int f;
	int sgn = 0;

	// HANDLE NAN AND INIFITY
	sgn = 0;
	if (isinf(ZREAL_(amax))) {
		// CHECK FOR NAN
		for (f = 0; f < N; f++) {
			if (isnan(ZREAL_(v[f]))) {
				ZSET_REAL_(amax, NAN);
				break;
			}
			if (isinf(ZREAL_(v[f]))) {
				// POSITIVE INFINITY
				if (signbit(ZREAL_(v[f])) == 0) {
					if (sgn < 0) { // INF - INF = NAN
						ZSET_REAL_(amax, NAN);
						break;
					}
					sgn = 1;
				}
				else {
					if (sgn > 0) { // INF - INF = NAN
						ZSET_REAL_(amax, NAN);
						break;
					}
					sgn = -1;
					ZSET_REAL_(amax, -INFINITY);
				}
			}
		}
	}
	sgn = 0;
	if (isinf(ZIMAG_(amax))) {
		// CHECK FOR NAN
		for (f = 0; f < N; f++) {
			if (isnan(ZIMAG_(v[f]))) {
				ZSET_IMAG_(amax, NAN);
				break;
			}
			if (isinf(ZIMAG_(v[f]))) {
				// POSITIVE INFINITY
				if (signbit(ZIMAG_(v[f])) == 0) {
					if (sgn < 0) { // INF - INF = NAN
						ZSET_IMAG_(amax, NAN);
						break;
					}
					sgn = 1;
				}
				else {
					if (sgn > 0) { // INF - INF = NAN
						ZSET_IMAG_(amax, NAN);
						break;
					}
					sgn = -1;
					ZSET_IMAG_(amax, -INFINITY);
				}
			}
		}
	}

	// NAN
	if (isnan(ZREAL_(amax))) {
		for (f = 0; f < fold; f++) 
			ZSET_REAL_(sum[f], ZREAL_(amax));
	}
	if (isnan(ZIMAG_(amax))) {
		for (f = 0; f < fold; f++) 
			ZSET_IMAG_(sum[f], ZIMAG_(amax));
	}

	// INFINITY
	if (isinf(ZREAL_(amax))) {
		if (!isinf(ZREAL_(sum[0]))) {
			for (f = 0; f < fold; f++) 
				ZSET_REAL_(sum[f], ZREAL_(amax));
		}
		if (signbit(ZREAL_(amax)) != signbit(ZREAL_(sum[0]))) {
			// INF - INF = NAN
			for (f = 0; f < fold; f++) 
				ZSET_REAL_(sum[f], NAN);
		}
	}
	// INFINITY
	if (isinf(ZIMAG_(amax))) {
		if (!isinf(ZIMAG_(sum[0]))) {
			for (f = 0; f < fold; f++) 
				ZSET_IMAG_(sum[f], ZIMAG_(amax));
		}
		if (signbit(ZIMAG_(amax)) != signbit(ZIMAG_(sum[0]))) {
			// INF - INF = NAN
			for (f = 0; f < fold; f++) 
				ZSET_IMAG_(sum[f], NAN);
		}
	}

	if (isnan(ZREAL_(sum[0])) && isnan(ZIMAG_(sum[0])))
		return -1;

	if (isinf(ZREAL_(sum[0])) && isinf(ZIMAG_(sum[0])))
		return 1;

	if (isnan(ZREAL_(sum[0])) && isinf(ZIMAG_(sum[0])))
		return 1;

	if (isinf(ZREAL_(sum[0])) && isnan(ZIMAG_(sum[0])))
		return 1;

	if (isnan(ZREAL_(sum[0])))
		return 2;

	if (isnan(ZIMAG_(sum[0])))
		return 3;

	return 0;
}

void zsumI1_(int N, int NB,
	double complex* v, int inc,
	int fold, double complex* sum, double complex* c) {

	double complex amax;

	int i, j;
	int status;
	int lN;
	int accu = 0;
	int maxN = dICapacity();
		
	for (i = 0; i < N; i+=NB, v+=NB*inc) {
		lN = NB < (N - i) ? NB:(N-i);

		// LOCAL MAX ABSOLUTE
		amax = zamax(lN, v, inc);
		status = zsum_exception(amax, N, v, inc, fold, sum);

		// TODO: CHECK STATUS
		if (status < 0)
			return;
		if (status == 1)
			continue;
		
		if (accu + lN > maxN) {
//			zIRenorm1(fold, sum, c, 1);
			accu = 0;
		}

		zmzupdate(&amax, sum, 1, c, 1, fold);

		zsumI2(lN, v, inc, fold, sum);

		accu += lN;

		if (status == 2) {
			for (j = 0; j < fold; j++)
				ZSET_REAL_(sum[j], ZREAL_(amax));
		}
		if (status == 3) {
			for (j = 0; j < fold; j++)
				ZSET_IMAG_(sum[j], ZIMAG_(amax));
		}
	}
	zIRenorm1(fold, sum, c, 1);
}

