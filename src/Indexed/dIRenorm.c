/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "Indexed.h"
#include "../Common/Common.h"

#ifdef __SSE2__
#	include <emmintrin.h>
#endif

void dIRenorm1(int n, double* X, double* leading, int inc) {
	int i;
	double M;
	double x;
	for (i = 0; i < n; i++, X += inc, leading += inc) {
		x = X[0];
		if (x == 0.0)
			continue;

		M = ufp(x);
		if (x >= (M * 1.75)) {
			X[0] -= M * 0.25;
			leading[0] += 1;
		}
		else if (x < (M * 1.25)) {
			X[0] += M * 0.5;
			leading[0] -= 2;
		}
		else if (x < (M * 1.5)) {
			X[0] += M * 0.25;
			leading[0] -= 1;
		}
	}
}

void dIRenorm2(int n, double* X, int incX) {
	dIRenorm1(n, X, X + n * incX, incX);
}

void zIRenorm1(int fold, double complex* rep, double complex* carry, int inc) {
	dIRenorm1(fold, (double*)rep, (double*)carry,2 * inc);
	dIRenorm1(fold, ((double*)rep) + 1, (double*)carry + 1, 2 * inc);
}


