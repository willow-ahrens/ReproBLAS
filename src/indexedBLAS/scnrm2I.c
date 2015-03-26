/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexedBLAS.h"
#include "../types.h"
#include "../config.h"

int rscnrm2_exception(
	float amax,
	int N, float complex* v, int incv,
	int fold, float complex* sum
) {
	int f;
		
	// NAN
	if (isnan(amax)) {
		for (f = 0; f < fold; f++) 
			CSET_(sum[f], amax, amax);
		return -1;
	}

	if (isinf(CREAL_(sum[0])) || isinf(CIMAG_(sum[0]))) {
		return 1;
	}

	// INFINITY
	if (isinf(amax)) {
		for (f = 0; f < fold; f++) 
			CSET_(sum[f], INFINITY, INFINITY);
		return 1;
	}

	return 0;
}

float scnrm2I1_(int N, int NB,
	float complex* v, int inc,
	int fold, int W, float* sum, F_CARRY_T* carry,
	float complex* work) {
	
	float amax;
	float complex amaxz;

	/*
	float complex BUFFER[MAX_FOLD];
	F_CARRY_T C[2*MAX_FOLD];
	*/

	int i, j;
	int status;
	int lN;
	int lB;
	int maxN = sICapacity();
	float scale = 0.0;
	float nscale = 0.0;
	int created = 0;

	if (work == NULL) {
		created = 1;
		work    = (float complex*) malloc(cISize(fold));
	}


	float complex* BUFFER = work;
	F_CARRY_T* C     = (F_CARRY_T*) (BUFFER + fold);

	for (i = 0; i < fold; i++) {
		CSET_(BUFFER[i], sum[i], 1.5 * ufpf(sum[i]));
	}
	for (i = 0; i < fold; i++) {
		C[2*i]= carry[i];
		C[2*i+1] = 0;
	}
		
	for (i = 0; i < N; i+=NB, v+=NB*inc) {
		lN = NB < (N - i) ? NB:(N-i);

		// LOCAL MAX ABSOLUTE
		amaxz = camax(lN, v, inc);
		amax = CREAL_(amaxz) < CIMAG_(amaxz) ? CIMAG_(amaxz) : CREAL_(amaxz);

		if (amax == 0.0)
			continue;

		status = rscnrm2_exception(amax, lN, v, inc, fold, BUFFER);

		if (status < 0) {
			for (j = 0; j < 2 * fold; j++)
				sum[i] = NAN;
			if (created == 1) {
				free(work);
			}
			return 1.0;
		}

		if (status > 0)
			continue;
		
		sIBoundary_(1, W, amax, &nscale, 1); 
		nscale = nscale / 1.5;
//		printf("max: %g, scale: %g, nscale: %g \n", amax, scale, nscale);

		// UPDATE NEW SCALE
		if (nscale > scale) {
			if (scale > 0.0) {
				scale = scale / nscale;
				scale = scale * scale;
				for (j = 0; j < fold; j++) {
					CSCAL_(BUFFER[j], BUFFER[j], scale);
				}
			}
			scale = nscale;

		}
			nscale = 1.0 / scale;
			amax *= nscale;
			amax = amax * amax;
//			printf("\nBefore Update: \n");
//			cIprint(fold, BUFFER);
			cIUpdates1(fold, W, BUFFER, C, 1, amax);

//			printf("\nAfter Update: \n");
//			cIprint(fold, BUFFER);


		for (j = 0; j < lN; j+=maxN) {
			lB = maxN < lN - j ? maxN : lN - j;
			scnrm2I2(lB, v+j*inc, inc, nscale, fold, BUFFER);

//			printf("\nBefore renorm: \n");
//			cIprint(fold, BUFFER);
			cIRenorm1(fold, BUFFER, C, 1);
//			printf("\nAfter renorm: \n");
//			cIprint(fold, BUFFER);
		}
//		printf("\nAfter 1 block: \n");
//		cIprint(fold, BUFFER);
		
	}

	for (i = 0; i < fold; i++) {
		sum[i] = (CREAL_(BUFFER[i]) - 1.5 * ufpf(CREAL_(BUFFER[i]))) + CIMAG_(BUFFER[i]);
	}
	for (i=0 ; i < fold; i++) {
		carry[i] = C[2*i]+ C[2*i+1];
	}

	sIRenorm1(fold, sum, carry, 1);

	if (created == 1) {
		free(work);
	}

	return scale;
}

