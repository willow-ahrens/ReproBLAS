/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

double dznrm2I(
	int N,
	double complex* x, int incx,
	I_double* sum
) {
	return dznrm2I1(N, x, incx, DEFAULT_FOLD, sum->m, sum->c);
}

double rdznrm2(int N, double complex* v, int inc) {
	I_double sum;
	dISetZero(sum);
	double scale = 0.0;
	scale = dznrm2I1(N, v, inc, DEFAULT_FOLD, sum.m, sum.c);

	if (isnan(scale))
		return scale;

	if (isinf(sum.m[0]))
		return sum.m[0];

	double sqrt_sum = ddiconv(&sum, DEFAULT_FOLD);
	sqrt_sum = sqrt(sqrt_sum);
	return scale * sqrt_sum;
}
