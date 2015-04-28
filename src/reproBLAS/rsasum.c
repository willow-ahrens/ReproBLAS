/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

I_float sasumI(int N, float* v, int inc) {
	I_float sum;

	sisetzero(DEFAULT_FOLD, &sum);
	sasumI1(N, v, inc, DEFAULT_FOLD, sum.m, sum.c);

	return sum;
}

float rsasum(int N, float* v, int inc) {
	I_float sum;

	sisetzero(DEFAULT_FOLD, &sum);
	sasumI1(N, v, inc, DEFAULT_FOLD, sum.m, sum.c);

	return ssiconv(&sum, DEFAULT_FOLD);
}

