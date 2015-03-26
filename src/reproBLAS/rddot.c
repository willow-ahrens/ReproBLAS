/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"



#define MACHEPS DBL_EPSILON
#define PREC    53
#define OVERFLOW_THRES DBL_MAX_EXP
#define CHECK_NAN_INF

I_double ddotI(
	int N,
	double* x, int incx,
	double* y, int incy
) {

	I_double dot;

	dISetZero(dot);

	ddotI1(N, x, incx, y, incy, DEFAULT_FOLD, 0, dot.m, dot.c);

	return dot;
}

double rddot(
	int N,
	double* x, int incx,
	double* y, int incy
) {

	I_double dot;

	dISetZero(dot);

	ddotI1(N, x, incx, y, incy, DEFAULT_FOLD, 0, dot.m, dot.c);

	return Iconv2d(dot);
}
