/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexedBLAS.h"

int rsnrm2_exception(
	float amax,
	int N, float* v, int incv,
	int fold, float* sum
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

float snrm2I1_(int N, int NB,
	float* v, int inc,
	int fold, int W, float* sum, float* c) {
	float amax;

	int i, j;
	int status;
	int lN;
	int maxN = sICapacity();
	float scale = 0.0;
	float nscale = 0.0;
		
	for (i = 0; i < N; i+=NB, v+=NB * inc) {
		lN = NB < (N - i) ? NB:(N-i);

		// LOCAL MAX ABSOLUTE
		amax = samax(lN, v, inc);

		if (amax == 0.0)
			continue;
#		ifdef CHECK_NAN_INF
		status = rsnrm2_exception(amax, lN, v, inc, fold, sum);
		if (status < 0)
			return 1.0;
		if (status > 0)
			continue;
#		endif

		sIBoundary_(1, W, amax, &nscale, 1); 
		nscale = ufpf(nscale);
//		nscale = nscale / 1.5;
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

		sIUpdate1(fold, W, amax, sum, c, 1);
		
		for (j = 0; j < lN - maxN + 1; j+=maxN) {
			snrm2I2(maxN, v + j * inc, inc, nscale, fold, sum);
			sIRenorm1(fold, sum, c, 1);
		}

		if (j < lN) {
			snrm2I2(lN - j, v + j * inc, inc, nscale, fold, sum);
			sIRenorm1(fold, sum, c, 1);
		}
	}

	return scale;
}

