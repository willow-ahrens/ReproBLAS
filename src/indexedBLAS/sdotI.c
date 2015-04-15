/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexedBLAS.h"

// TODO CHECK FOR FALSE INFINITY WHEN ABSOLUTE MAX IS CLOSE TO INFINITY
// SOLUTION: SCALE INPUTS

int rsdot_exception(
	float amax,
	int N, 
	float* x, int incx,
	float* y, int incy,
	int fold, float* dot
) {
	float dindex;

	int f;
	int sgn = 0;


	// HANDLE NAN AND INIFITY
	sgn = 0;
	if (isinf(amax)) {
		// CHECK FOR NAN
		float ptmp;
		for (f = 0; f < N; f++) {
			ptmp = x[f] * y[f];
			if (isnan(ptmp)) {
				amax = NAN;
				break;
			}
			if (ptmp) {
				// POSITIVE INFINITY
				if (signbit(ptmp) == 0) {
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

	if (isnan(amax)) {
		for (f = 0; f < fold; f++) 
			dot[f] = amax;
#ifdef DEBUG
		fprintf(stdout, "EARLY TERMINATION DUE TO NAN EXECPTION AT INDEX %d \n", f);
#endif
		return -1;
	}

	if (isinf(amax)) {
		if (!isinf(dot[0])) {
			for (f = 0; f < fold; f++) 
				dot[f] = amax;
			return 1;
		}
		if (signbit(amax) == signbit(dot[1])) {
			return 1;
		}
		// INF - INF = NAN
		for (f = 0; f < fold; f++) 
			dot[f] = NAN;
#ifdef DEBUG
		fprintf(stdout, "EARLY TERMINATION DUE TO NAN EXECPTION AT %d \n", f);
#endif
		return -1;
	}

	return 0;
}

void sdotI1_(
	int N, int NB,
	float* x, int incx,
	float* y, int incy,
	int fold, 
	float* dot, float* c
) {

	int lN, i, j;
	int status;
	float amax;
	int maxN = sICapacity();

	for (i = 0; i < N; i+=NB, x += NB*incx, y += NB*incy) {
		lN = NB < (N - i) ? NB : N - i;
		amax = samaxm(lN, x, incx, y, incy);

		if (amax == 0.0)
			continue;
#		ifdef CHECK_NAN_INF
		status = rsdot_exception(amax, lN, x, incx, y, incy, fold, dot);
		if (status < 0)
			return;
		if (status > 0)
			continue;
#		endif

		sIUpdate1(fold, amax, dot, c, 1);
		
		for (j = 0; j < lN - maxN + 1; j+=maxN) {
			sdotI2(maxN, x + j * incx, incx, y + j * incy, incy, fold, dot);
			sIRenorm1(fold, dot, c, 1);
		}

		if (j < lN) {
			sdotI2(lN - j, x + j * incx, incx, y + j * incy, incy, fold, dot);
			sIRenorm1(fold, dot, c, 1);
		}
	}
}

