/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

I_float_Complex csumI(int N, float complex* v, int inc) {
	I_float_Complex sum;
	cisetzero(DEFAULT_FOLD, &sum);
	csumI1(N, v, inc, DEFAULT_FOLD, (float complex*)sum.m, sum.c);
	return sum;
}

float complex rcsum(int N, float complex* v, int inc) {
	I_float_Complex sum;
	cisetzero(DEFAULT_FOLD, &sum);
    float complex ret;
	csumI1(N, v, inc, DEFAULT_FOLD, (float complex*)sum.m, sum.c);
    cciconv_sub(&sum, &ret, DEFAULT_FOLD);
    return ret;
}

