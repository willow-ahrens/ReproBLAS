/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rblas1.h"

// TODO CHECK FOR FALSE INFINITY WHEN ABSOLUTE MAX IS CLOSE TO INFINITY
// SOLUTION: SCALE INPUTS

int exceptionHandling(
	double amax,
	int N, double* v, int incv,
	int fold, double* sum
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

void dsumI1_(int N, int NB,
	double* v, int inc,
	int fold, int W,
	double* sum, double* c) {
	
	double amax;
	int i;
	int status;
	int lN;
	int accu = 0;
	int maxN = dICapacity();
		
	for (i = 0; i < N; i+=NB, v+=NB) {
		lN = NB < (N - i) ? NB:(N-i);

		// LOCAL MAX ABSOLUTE
		amax = DAMAX(lN, v, inc);
		if (amax == 0.0)
			continue;

		status = exceptionHandling(amax, N, v, inc, fold, sum);

		if (status < 0)
			break;
		if (status > 0)
			continue;

		if (accu + lN > maxN) {
			dIRenorm1(fold, sum, c, 1);
			accu = 0;
		}
		
		dIUpdate1(fold, W, sum, c, 1, amax);
		
		// TODO: CHECK POTENTIAL FALSE INFINITY
		dsumI2(lN, v, inc, fold, sum);
		accu += lN;
	}
	dIRenorm1(fold, sum, c, 1);
}

