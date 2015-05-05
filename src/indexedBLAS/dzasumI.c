/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexedBLAS.h"
#include "../config.h"
#include "../types.h"

int rdzasum_exception(
	double amax,
	int N, dcomplex* v, int incv,
	int fold, dcomplex* sum
) {
	int f;
		
	// NAN
	if (isnan(amax)) {
		for (f = 0; f < fold; f++) 
			ZSET_(sum[f], amax, amax);
#ifdef DEBUG
		fprintf(stdout, "EARLY TERMINATION DUE TO NAN EXECPTION AT INDEX %d \n", f);
#endif
		return -1;
	}

	if (isinf(ZREAL_(sum[0])) || isinf(ZIMAG_(sum[0]))) {
		return 1;
	}

	// INFINITY
	if (isinf(amax)) {
		for (f = 0; f < fold; f++) 
			ZSET_(sum[f], INFINITY, INFINITY);
		return 1;
	}

	return 0;
}

void dzasumI1_(int N, int NB,
	dcomplex* v, int inc,
	int fold, double* sum, double * c) {

    dmzasum(fold, N, v, inc, sum, 1, c, 1);
    return;

	double amax;
	double M;
	dcomplex amaxz;

	int i, j;
	int status;
	int lN;
	int accu = 0;
	int maxN = dicapacity();
	double tmp;
	dcomplex BUFFER[2 * MAX_FOLD];

	for (i = 0; i < fold; i++) {
		ZSET_(BUFFER[i], sum[i], 1.5*ufp(sum[i]));
		ZSET_(BUFFER[i+fold], c[i], 0.0);
	}
		
	for (i = 0; i < N; i+=NB, v+=NB*inc) {
		lN = NB < (N - i) ? NB:(N-i);

		// LOCAL MAX ABSOLUTE
		zamax_sub(lN, v, inc, &amaxz);
		amax = ZREAL_(amaxz) < ZIMAG_(amaxz) ? ZIMAG_(amaxz) : ZREAL_(amaxz);

		if (amax == 0.0)
			continue;

		status = rdzasum_exception(amax, lN, v, inc, fold, BUFFER);

//		printf("max: %g, status: %d\n", amax, status);

		if (status < 0) {
			for (j = 0; j < 2 * fold; j++)
				sum[j] = NAN;
			return;
		}

		if (status > 0)
			continue;
		
		zmdupdate(fold, amax, BUFFER, 1, BUFFER + fold, 1);

		if (accu + lN > maxN) {
			zmrenorm(fold, BUFFER, 1, BUFFER + fold, 1);
			accu = 0;
		}

		// TODO: CHECK POTENTIAL FALSE INFINITY
		dzasumI2(lN, v, inc, fold, BUFFER);

		accu += lN;
	}

	zmrenorm(fold, BUFFER, 1, BUFFER + fold, 1);

	for (i = 0; i < fold; i++) {
		sum[i] = (ZREAL_(BUFFER[i])- 1.5*ufp(ZIMAG_(BUFFER[i]))) + (ZIMAG_(BUFFER[i]) );
	}
	for (; i < 2 * fold; i++) {
		sum[i] = ZREAL_(BUFFER[i]) + ZIMAG_(BUFFER[i]);
	}
	dmrenorm(fold, sum, 1, c, 1);
}

