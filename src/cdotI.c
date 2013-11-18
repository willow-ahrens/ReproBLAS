/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rblas1.h"
#include "config.h"
#include "types.h"

// TODO CHECK FOR FALSE INFINITY WHEN ABSOLUTE MAX IS CLOSE TO INFINITY
// SOLUTION: SCALE INPUTS

int cdot_exception(
	float complex amax,
	int N,
	float complex* v, int incv,
	float complex* y, int incy,
	int fold, float complex* dot
) {
	int f;
	int sgn = 0;

	// HANDLE NAN AND INIFITY
	sgn = 0;
	if (isinf(CREAL_(amax))) {
		// CHECK FOR NAN
		for (f = 0; f < N; f++) {
			if (isnan(CREAL_(v[f]))) {
				CSET_REAL_(amax, NAN);
				break;
			}
			if (isinf(CREAL_(v[f]))) {
				// POSITIVE INFINITY
				if (signbit(CREAL_(v[f])) == 0) {
					if (sgn < 0) { // INF - INF = NAN
						CSET_REAL_(amax, NAN);
						break;
					}
					sgn = 1;
				}
				else {
					if (sgn > 0) { // INF - INF = NAN
						CSET_REAL_(amax, NAN);
						break;
					}
					sgn = -1;
					CSET_REAL_(amax, -INFINITY);
				}
			}
		}
	}
	sgn = 0;
	if (isinf(CIMAG_(amax))) {
		// CHECK FOR NAN
		for (f = 0; f < N; f++) {
			if (isnan(CIMAG_(v[f]))) {
				CSET_IMAG_(amax, NAN);
				break;
			}
			if (isinf(CIMAG_(v[f]))) {
				// POSITIVE INFINITY
				if (signbit(CIMAG_(v[f])) == 0) {
					if (sgn < 0) { // INF - INF = NAN
						CSET_IMAG_(amax, NAN);
						break;
					}
					sgn = 1;
				}
				else {
					if (sgn > 0) { // INF - INF = NAN
						CSET_IMAG_(amax, NAN);
						break;
					}
					sgn = -1;
					CSET_IMAG_(amax, -INFINITY);
				}
			}
		}
	}

	// NAN
	if (isnan(CREAL_(amax))) {
		for (f = 0; f < fold; f++) 
			CSET_REAL_(dot[f], CREAL_(amax));
	}
	if (isnan(CIMAG_(amax))) {
		for (f = 0; f < fold; f++) 
			CSET_IMAG_(dot[f], CIMAG_(amax));
	}

	// INFINITY
	if (isinf(CREAL_(amax))) {
		if (!isinf(CREAL_(dot[0]))) {
			for (f = 0; f < fold; f++) 
				CSET_REAL_(dot[f], CREAL_(amax));
		}
		if (signbit(CREAL_(amax)) != signbit(CREAL_(dot[0]))) {
			// INF - INF = NAN
			for (f = 0; f < fold; f++) 
				CSET_REAL_(dot[f], NAN);
		}
	}
	// INFINITY
	if (isinf(CIMAG_(amax))) {
		if (!isinf(CIMAG_(dot[0]))) {
			for (f = 0; f < fold; f++) 
				CSET_IMAG_(dot[f], CIMAG_(amax));
		}
		if (signbit(CIMAG_(amax)) != signbit(CIMAG_(dot[0]))) {
			// INF - INF = NAN
			for (f = 0; f < fold; f++) 
				CSET_IMAG_(dot[f], NAN);
		}
	}

	if (isnan(CREAL_(dot[0])) && isnan(CIMAG_(dot[0])))
		return -1;

	if (isinf(CREAL_(dot[0])) && isinf(CIMAG_(dot[0])))
		return 1;

	if (isnan(CREAL_(dot[0])) && isinf(CIMAG_(dot[0])))
		return 1;

	if (isinf(CREAL_(dot[0])) && isnan(CIMAG_(dot[0])))
		return 1;

	if (isnan(CREAL_(dot[0])))
		return 2;

	if (isnan(CIMAG_(dot[0])))
		return 3;

	return 0;
}

void cdotI1_(int N, int NB,
	float complex* v, int inc,
	float complex* y, int incy,
	int fold, int W, float complex* dot, F_CARRY_T* C,
	int conj) {
/*
 * Purpose
 *  Compute the inner product (dot product)
 *  two single complex vectors
 * 
 * Arguments
 *  N     [input]  size of input vectors
 *  NB    [input]  Block size of input vectors
 *  v     [input]  pointer of the first element of the first array
 *  inc   [input]  increment between two consecutive elements of V
 *  y     [input]  pointer of the first element of the second array
 *  incx  [input]  increment between two consecutive elements of Y
 *  fold  [input]  Number of bins of Indexed Format
 *  W     [input]  Bin width of Indexed Format
 *  dot   [in/out] Mantissa array in Indexed Format of the dot product
 *  C     [in/out] Carry-bits array in Indexed Format of the dot product
 *  conj  [input]  whether Y is conjugated or not
 *                 0: not conjugated    1: conjugated
 */

	float complex amax;

	int i, j;
	int status;
	int lN;
	int maxN = sICapacity() / 2;
	int lB;
		
	for (i = 0; i < N; i+=NB, v+=NB*inc, y+=NB*incy) {
		lN = NB < (N - i) ? NB:(N-i);

		// LOCAL MAX ABSOLUTE
		amax = camaxm(lN, v, inc, y, incy);
//		printf("\ndotI1, max: 2^%g 2^%g\n", log2f(amax.real), log2f(amax.imag)); //		amax = camax(lN, v, inc);
		status = cdot_exception(amax, N, v, inc, y, incy, fold, dot);

		// TODO: CHECK STATUS
		if (status < 0)
			return;
		if (status == 1)
			continue;
		
		cIUpdate1(fold, W, dot, C, 1, amax);

		for (j = 0; j < lN; j+=maxN) {
			lB = maxN < lN - j ? maxN : lN - j;
			if (conj)
				cdotcI2(lB, v+j*inc, inc, y+j*incy, incy, fold, dot);
			else 
				cdotuI2(lB, v+j*inc, inc, y+j*incy, incy, fold, dot);

			cIRenorm1(fold, dot, C, 1);
		}

		if (status == 2) {
			for (j = 0; j < fold; j++)
				CSET_REAL_(dot[j], CREAL_(amax));
		}
		if (status == 3) {
			for (j = 0; j < fold; j++)
				CSET_IMAG_(dot[j], CIMAG_(amax));
		}

	}
}

