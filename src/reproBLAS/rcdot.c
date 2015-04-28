/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

I_float_Complex cdotcI(
	int N,
	float complex* x, int incx,
	float complex* y, int incy
) {
	I_float_Complex dot;
	cisetzero(DEFAULT_FOLD, &dot);
	cdotcI1(N, x, incx, y, incy, DEFAULT_FOLD, (float complex*)dot.m, dot.c);
	return dot;
}

I_float_Complex cdotuI(
	int N,
	float complex* x, int incx,
	float complex* y, int incy
) {
	I_float_Complex dot;
	cisetzero(DEFAULT_FOLD, &dot);
	cdotuI1(N, x, incx, y, incy, DEFAULT_FOLD, (float complex*)dot.m, dot.c);
	return (dot);
}

float complex rcdotc(
	int N,
	float complex* x, int incx,
	float complex* y, int incy
) {
	I_float_Complex dot;
	cisetzero(DEFAULT_FOLD, &dot);
    float complex ret;
	cdotcI1(N, x, incx, y, incy, DEFAULT_FOLD, (float complex*)dot.m, dot.c);
    cciconv_sub(DEFAULT_FOLD, &dot, &ret);
    return ret;
}

float complex rcdotu(
	int N,
	float complex* x, int incx,
	float complex* y, int incy
) {
	I_float_Complex dot;
	cisetzero(DEFAULT_FOLD, &dot);
    float complex ret;
	cdotuI1(N, x, incx, y, incy, DEFAULT_FOLD, (float complex*)dot.m, dot.c);
    cciconv_sub(DEFAULT_FOLD, &dot, &ret);
    return ret;
}

