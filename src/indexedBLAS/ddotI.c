/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexedBLAS.h"

int rddot_exception(
	double amax,
	int N, 
	double* x, int incx,
	double* y, int incy,
	int fold, double* dot
) {
	double dindex;

	int f;
	int sgn = 0;


	// HANDLE NAN AND INIFITY
	sgn = 0;
	if (isinf(amax)) {
		// CHECK FOR NAN
		double ptmp;
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

void ddotI1_(
	int N, int NB,
	double* x, int incx,
	double* y, int incy,
	int fold, double* dot, double* c
) {

	int lN, i;
	int status;
	double amax;
	int accu = 0;
	int maxN = dICapacity();

	for (i = 0; i < N; i+=NB, x += NB*incx, y += NB*incy) {
		// BLOCK SIZE
		lN = NB < (N - i) ? NB : N - i;

		// LOCAL MAXIMUM ABSOLUTE VALUE
		amax = damaxm(lN, x, incx, y, incy);

		if (amax == 0.0)
			continue;
#		ifdef CHECK_NAN_INF
		status = rddot_exception(amax, lN, x, incx, y, incy, fold, dot);
		if (status < 0)
			return;
		if (status > 0)
			continue;
#		endif

		if (accu + lN > maxN) {
			dmrenorm(dot, 1, c, 1, fold);
			accu = 0;
		}
		
		dmdupdate(amax, dot, 1, c, 1, fold);
		
		ddotI2(lN, x, incx, y, incy, fold, dot);

		accu += lN;
	}
	dmrenorm(dot, 1, c, 1, fold);
}

