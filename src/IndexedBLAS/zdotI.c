/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "IndexedBLAS.h"
#include "../types.h"

// TODO CHECK FOR FALSE INFINITY WHEN ABSOLUTE MAX IS CLOSE TO INFINITY
// SOLUTION: SCALE INPUTS

int zdot_exception(
	double complex amax,
	int N,
	double complex* v, int incv,
	double complex* y, int incy,
	int fold, double complex* dot
) {
	int f;
	int sgn = 0;

	// HANDLE NAN AND INIFITY
	sgn = 0;
	if (isinf(ZREAL_(amax))) {
		// CHECK FOR NAN
		for (f = 0; f < N; f++) {
			if (isnan(ZREAL_(v[f]))) {
				ZSET_REAL_(amax, NAN);
				break;
			}
			if (isinf(ZREAL_(v[f]))) {
				// POSITIVE INFINITY
				if (signbit(ZREAL_(v[f])) == 0) {
					if (sgn < 0) { // INF - INF = NAN
						ZSET_REAL_(amax, NAN);
						break;
					}
					sgn = 1;
				}
				else {
					if (sgn > 0) { // INF - INF = NAN
						ZSET_REAL_(amax, NAN);
						break;
					}
					sgn = -1;
					ZSET_REAL_(amax, -INFINITY);
				}
			}
		}
	}
	sgn = 0;
	if (isinf(ZIMAG_(amax))) {
		// CHECK FOR NAN
		for (f = 0; f < N; f++) {
			if (isnan(ZIMAG_(v[f]))) {
				ZSET_IMAG_(amax, NAN);
				break;
			}
			if (isinf(ZIMAG_(v[f]))) {
				// POSITIVE INFINITY
				if (signbit(ZIMAG_(v[f])) == 0) {
					if (sgn < 0) { // INF - INF = NAN
						ZSET_IMAG_(amax, NAN);
						break;
					}
					sgn = 1;
				}
				else {
					if (sgn > 0) { // INF - INF = NAN
						ZSET_IMAG_(amax, NAN);
						break;
					}
					sgn = -1;
					ZSET_IMAG_(amax, -INFINITY);
				}
			}
		}
	}

	// NAN
	if (isnan(ZREAL_(amax))) {
		for (f = 0; f < fold; f++) 
			ZSET_REAL_(dot[f], ZREAL_(amax));
	}
	if (isnan(ZIMAG_(amax))) {
		for (f = 0; f < fold; f++) 
			ZSET_IMAG_(dot[f], ZIMAG_(amax));
	}

	// INFINITY
	if (isinf(ZREAL_(amax))) {
		if (!isinf(ZREAL_(dot[0]))) {
			for (f = 0; f < fold; f++) 
				ZSET_REAL_(dot[f], ZREAL_(amax));
		}
		if (signbit(ZREAL_(amax)) != signbit(ZREAL_(dot[0]))) {
			// INF - INF = NAN
			for (f = 0; f < fold; f++) 
				ZSET_REAL_(dot[f], NAN);
		}
	}
	// INFINITY
	if (isinf(ZIMAG_(amax))) {
		if (!isinf(ZIMAG_(dot[0]))) {
			for (f = 0; f < fold; f++) 
				ZSET_IMAG_(dot[f], ZIMAG_(amax));
		}
		if (signbit(ZIMAG_(amax)) != signbit(ZIMAG_(dot[0]))) {
			// INF - INF = NAN
			for (f = 0; f < fold; f++) 
				ZSET_IMAG_(dot[f], NAN);
		}
	}

	if (isnan(ZREAL_(dot[0])) && isnan(ZIMAG_(dot[0])))
		return -1;

	if (isinf(ZREAL_(dot[0])) && isinf(ZIMAG_(dot[0])))
		return 1;

	if (isnan(ZREAL_(dot[0])) && isinf(ZIMAG_(dot[0])))
		return 1;

	if (isinf(ZREAL_(dot[0])) && isnan(ZIMAG_(dot[0])))
		return 1;

	if (isnan(ZREAL_(dot[0])))
		return 2;

	if (isnan(ZIMAG_(dot[0])))
		return 3;

	return 0;
}

void zdotI1_(int N, int NB,
	double complex* v, int inc,
	double complex* y, int incy,
	int fold, int W, double complex* dot, double complex* c,
	int conj) {

	double complex amax;

	int i, j;
	int status;
	int lN;
	int accu = 0;
	int maxN = dICapacity();
		
	for (i = 0; i < N; i+=NB, v+=NB*inc, y+=NB*incy) {
    double complex foo = *y;
		lN = NB < (N - i) ? NB:(N-i);

		// LOCAL MAX ABSOLUTE
		amax = zamaxm(lN, v, inc, y, incy);
		status = zdot_exception(amax, N, v, inc, y, incy, fold, dot);

		// TODO: CHECK STATUS
		if (status < 0)
			return;
		if (status == 1)
			continue;
		
		if (accu + lN > maxN) {
			zIRenorm1(fold, dot, c, 1);
			accu = 0;
		}

		zIUpdate1(fold, W, dot, c, 1, amax);

		if (conj)
			zdotcI2(lN, v, inc, y, incy, fold, dot);
		else 
			zdotuI2(lN, v, inc, y, incy, fold, dot);
			

		accu += lN;

		if (status == 2) {
			for (j = 0; j < fold; j++)
				ZSET_REAL_(dot[j], ZREAL_(amax));
		}
		if (status == 3) {
			for (j = 0; j < fold; j++)
				ZSET_IMAG_(dot[j], ZIMAG_(amax));
		}

	}
	zIRenorm1(fold, dot, c, 1);
}

