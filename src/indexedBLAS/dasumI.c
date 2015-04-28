/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexedBLAS.h"

int rdasum_exception(
	double amax,
	int N, double* v, int incv,
	int fold, double* sum
) {
	int f;
		
	// NAN
	if (isnan(amax)) {
		for (f = 0; f < fold; f++) 
			sum[f] = amax;
#ifdef DEBUG
		fprintf(stdout, "EARLY TERMINATION DUE TO NAN EXECPTION AT INDEX %d \n", f);
#endif
		return -1;
	}

	if (isinf(sum[0])) {
		return 1;
	}

	// INFINITY
	if (isinf(amax)) {
		for (f = 0; f < fold; f++) 
			sum[f] = INFINITY;
		return 1;
	}

	return 0;
}

void dasumI1_(int N, int NB,
	double* v, int inc,
	int fold, 
	double* sum, double* c) {

	double amax;
	int i;
	int status;
	int lN;
	int accu = 0;
	int maxN = dicapacity();
		
	for (i = 0; i < N; i+=NB, v+=NB * inc) {
		lN = NB < (N - i) ? NB:(N-i);

		// LOCAL MAX ABSOLUTE
		amax = damax(lN, v, inc);
		status = rdasum_exception(amax, lN, v, inc, fold, sum);

//		printf("max: %g, status: %d\n", amax, status);

		if (status < 0)
			return;

		if (status > 0)
			continue;
		
		if (accu + lN > maxN) {
			dmrenorm(sum, 1, c, 1, fold);
			accu = 0;
		}

		dmdupdate(amax, sum, 1, c, 1, fold);
		
		// TODO: CHECK POTENTIAL FALSE INFINITY
		dasumI2(lN, v, inc, fold, sum);

		accu += lN;
	}
	dmrenorm(sum, 1, c, 1, fold);
}

