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

int rdznrm2_exception(
	double amax,
	int N, double complex* v, int incv,
	int fold, double complex* sum
) {
	int f;
		
	// NAN
	if (isnan(amax)) {
		for (f = 0; f < fold; f++) 
			ZSET_(sum[f], amax, amax);
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

double dznrm2I1_(int N, int NB,
	double complex* v, int inc,
	int fold, double* sum, double* c) {
	double amax;
	double complex amaxz;
	double complex BUFFER[2 * MAX_FOLD];

	int i, j;
	int status;
	int lN;
	int accu = 0;
	int maxN = dICapacity();
	double scale = 0.0;
	double nscale = 0.0;
	double M;

	for (i = 0; i < fold; i++) {
		ZSET_(BUFFER[i], sum[i], sum[i]);
		ZSET_(BUFFER[i+fold], c[i], 0);
	}
		
	for (i = 0; i < N; i+=NB, v+=NB*inc) {
		lN = NB < (N - i) ? NB:(N-i);

		// LOCAL MAX ABSOLUTE
		amaxz = zamax(lN, v, inc);
		amax = ZREAL_(amaxz) < ZIMAG_(amaxz) ? ZIMAG_(amaxz) : ZREAL_(amaxz);

		if (amax == 0.0)
			continue;

		status = rdznrm2_exception(amax, lN, v, inc, fold, BUFFER);

		if (status < 0) {
			for (j = 0; j < 2 * fold; j++)
				sum[i] = NAN;
			return 1.0;
		}

		if (status > 0)
			continue;
		
		if (accu + lN > maxN) {
			zIRenorm1(fold, BUFFER, BUFFER + fold, 1);
			accu = 0;
		}

		dIBoundary(1, amax, &nscale, 1); 
		nscale = nscale / 1.5;
//		printf("max: %g, nscale: %g \n", amax, nscale);

		// UPDATE NEW SCALE
		if (nscale > scale) {
			if (scale > 0.0) {
				scale = scale / nscale;
				scale = scale * scale;
				for (j = 0; j < fold; j++) {
					ZSCAL_(BUFFER[j], BUFFER[j], scale);
				}
			}
			scale = nscale;
		}


		nscale = 1.0 / scale;
		amax *= nscale;
		amax = amax * amax;

		zIUpdates1(fold, BUFFER, BUFFER + fold, 1, amax);
		
		// TODO: CHECK POTENTIAL FALSE INFINITY
		dznrm2I2(lN, v, inc, nscale, fold, BUFFER);

		accu += lN;
	}
	
	zIRenorm1(fold, BUFFER, BUFFER + fold, 1);

	for (i = 0; i < fold; i++) {
		sum[i] = (ZREAL_(BUFFER[i]) - sum[i]) + ZIMAG_(BUFFER[i]);
		c[i]   = ZREAL_(BUFFER[i+fold]) + ZIMAG_(BUFFER[i+fold]);
	}

	dIRenorm1(fold, sum, c,1);

	return scale;
}

