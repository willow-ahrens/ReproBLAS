/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

I_float scasumI(int N, float complex* v, int inc) {
	I_float sum;
	I_float_Complex tmp;
	sisetzero(DEFAULT_FOLD, &sum);
	scasumI1(N, v, inc, DEFAULT_FOLD, sum.m, sum.c, (float complex*)&tmp);
	return sum;
}

float rscasum(int N, float complex* v, int inc) {
	I_float sum;
	I_float_Complex tmp;
	sisetzero(DEFAULT_FOLD, &sum);
	scasumI1(N, v, inc, DEFAULT_FOLD, sum.m, sum.c, (float complex*)&tmp);
	return ssiconv(DEFAULT_FOLD, &sum);
}

