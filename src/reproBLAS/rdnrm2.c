/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

double dnrm2I(
	int N,
	double* x, int incx,
	I_double* sum
) {
	return dnrm2I1(N, x, incx, DEFAULT_FOLD, 0, sum->m, sum->c);
}

double rdnrm2(int N, double* v, int inc) {
	I_double sum;
	double sqrt_sum;
	double scale;

	dISetZero(sum);

	scale = dnrm2I1(N, v, inc, DEFAULT_FOLD, 0, sum.m, sum.c);

	if (isnan(scale))
		return scale;

	if (isinf(sum.m[0]))
		return sum.m[0];

	sqrt_sum = Iconv2d(sum);
	sqrt_sum = sqrt(sqrt_sum);
	return scale * sqrt_sum;
}
