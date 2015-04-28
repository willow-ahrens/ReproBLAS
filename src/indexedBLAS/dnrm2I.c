/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexedBLAS.h"

int rdnrm2_exception(
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

double dnrm2I1_(int N, int NB,
	double* v, int inc,
	int fold, double* sum, double* c
) {
	double amax;

	int i, j;
	int status;
	int lN;
	int accu = 0;
	int maxN = dicapacity();
	double scale = 0.0;
	double nscale = 0.0;
		
	for (i = 0; i < N; i+=NB, v+=NB * inc) {
		lN = NB < (N - i) ? NB:(N-i);

		// LOCAL MAX ABSOLUTE
		amax = damax(lN, v, inc);

		if (amax == 0.0)
			continue;

		status = rdnrm2_exception(amax, lN, v, inc, fold, sum);

		if (status < 0)
			return 1.0;

		if (status > 0)
			continue;
		
		if (accu + lN > maxN) {
			dmrenorm(sum, 1, c, 1, fold);
			accu = 0;
		}

		nscale = dbound(dindex(amax))/1.5; 
//		printf("max: %g, nscale: %g \n", amax, nscale);

		// UPDATE NEW SCALE
		if (nscale > scale) {
			if (scale > 0.0) {
				scale = scale / nscale;
				scale = scale * scale;
				for (j = 0; j < fold; j++) {
					sum[j] *= scale;
				}
			}
			scale = nscale;
		}


		nscale = 1.0 / scale;
		amax *= nscale;
		amax = amax * amax;

		dmdupdate(amax, sum, 1, c, 1, fold);
		
		// TODO: CHECK POTENTIAL FALSE INFINITY
		dnrm2I2(lN, v, inc, nscale, fold, sum);

		accu += lN;
	}
	dmrenorm(sum, 1, c, 1, fold);
	return scale;
}

