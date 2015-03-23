/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "IndexedBLAS.h"

int rssum_exception(
	float amax,
	int N, float* v, int incv,
	int fold, float* sum
) {
	int f;
	int sgn = 0;
		
	// HANDLE NAN AND INIFITY
	sgn = 0;
	if (isinf(amax)) {
		// CHECK FOR NAN
		for (f = 0; f < N; f++) {
			if (isnan(v[f])) {
				amax = NAN;
				break;
			}
			if (isinf(v[f])) {
				// POSITIVE INFINITY
				if (signbit(v[f]) == 0) {
					if (sgn < 0) { // INF - INF = NAN
						amax = NAN;
						break;
					}
					sgn = 1;
				}
				else {
					if (sgn > 0) { // INF - INF = NAN
						amax = NAN;
						break;
					}
					sgn = -1;
					amax = -INFINITY;
				}
			}
		}
	}

	// NAN
	if (isnan(amax)) {
		for (f = 0; f < fold; f++) 
			sum[f] = amax;
#ifdef DEBUG
		fprintf(stdout, "EARLY TERMINATION DUE TO NAN EXECPTION AT INDEX %d \n", i);
#endif
		return -1;
	}

	// INFINITY
	if (isinf(amax)) {
		if (!isinf(sum[0])) {
			for (f = 0; f < fold; f++) 
				sum[f] = amax;
			return 1;
		}
		if (signbit(amax) == signbit(sum[1])) {
			return 1;
		}
		// INF - INF = NAN
		for (f = 0; f < fold; f++) 
			sum[f] = NAN;
#ifdef DEBUG
		fprintf(stdout, "EARLY TERMINATION DUE TO NAN EXECPTION AT %d \n", i);
#endif
		return -1;
	}

	return 0;
}

void ssumI1_(int N, int NB,
	float* v, int inc,
	int fold, int W, float* sum, F_CARRY_T* c) {
	float amax;

	int i, j;
	int status;
	int lN;
	int maxN = sICapacity();
		
	for (i = 0; i < N; i+=NB, v+=NB * inc) {
		lN = NB < (N - i) ? NB:(N-i);

		// LOCAL MAX ABSOLUTE
		amax = samax(lN, v, inc);

		if (amax == 0.0)
			continue;
#		ifdef CHECK_NAN_INF
		status = rssum_exception(amax, lN, v, inc, fold, sum);
		if (status < 0)
			return;
		if (status > 0)
			continue;
#		endif

		sIUpdate1(fold, W, amax, sum, c, 1);
		
		for (j = 0; j < lN - maxN + 1; j+=maxN) {
			ssumI2(maxN, v + j*inc, inc, fold, sum);
			sIRenorm1(fold, sum, c, 1);
		}

		if (j < lN) {
			ssumI2(lN - j, v + j*inc, inc, fold, sum);
			sIRenorm1(fold, sum, c, 1);
		}
	}
}

