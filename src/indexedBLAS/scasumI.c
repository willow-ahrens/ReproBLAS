/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "../config.h"
#include "indexedBLAS.h"
#include "../types.h"

int rscasum_exception(
	float amax,
	int N, float complex* v, int incv,
	int fold, float complex* sum
) {
	int f;
		
	// NAN
	if (isnan(amax)) {
		for (f = 0; f < fold; f++) 
			CSET_(sum[f], amax, amax);
#ifdef DEBUG
		fprintf(stdout, "EARLY TERMINATION DUE TO NAN EXECPTION AT INDEX %d \n", f);
#endif
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

// work is a buffer of size fold (sizeof(float) + sizeof (carry))
void scasumI1_(int N, int NB,
	float complex* v, int inc,
	int fold, float* sum, float* carry,
	float complex* work) {

	float amax;
	float complex amaxz;

	int i, j;
	int status;
	int lN;
	int lB;
	int maxN = sicapacity();
	int created = 0;

	/*
	float complex BUFFER[MAX_FOLD];
	float C[2*MAX_FOLD];
	*/

	if (work == NULL) {
		created = 1;
		work    = (float complex*) malloc(cisize(fold));
	}

	float complex* BUFFER = work;
	float* C     = (float*) (BUFFER + fold);

	for (i = 0; i < fold; i++) {
		CSET_(BUFFER[i], sum[i], 1.5*ufpf(sum[i]));
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

		status = rscasum_exception(amax, lN, v, inc, fold, BUFFER);

//		printf("max: %g, status: %d\n", amax, status);

		if (status < 0) {
			for (j = 0; j < 2 * fold; j++)
				sum[j] = NAN;
			if (created == 1) {
				free(work);
			}
			return;
		}

		if (status > 0)
			continue;
		
		cmsupdate(fold, amax, BUFFER, 1, C, 1);

		for (j = 0; j < lN; j+=maxN) {
			lB = maxN < lN - j ? maxN : lN - j;
			scasumI2(lB, v+j*inc, inc, fold, BUFFER);

			cmrenorm(fold, BUFFER, 1, C, 1);
		}

	}

	for (i = 0; i < fold; i++) {
		sum[i] = (CREAL_(BUFFER[i]) - 1.5*ufpf(CREAL_(BUFFER[i]))) + CIMAG_(BUFFER[i]);
	}
	for (i=0 ; i < fold; i++) {
		carry[i] = C[2*i]+ C[2*i+1];
	}
	smrenorm(fold, sum, 1, carry, 1);

	if (created == 1) {
		free(work);
	}
}

