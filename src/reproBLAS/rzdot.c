#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

I_double_Complex zdotcI(
	int N,
	double complex* x, int incx,
	double complex* y, int incy
) {
	I_double_Complex dot;
	zisetzero(DEFAULT_FOLD, &dot);
	zdotcI1(N, x, incx, y, incy, DEFAULT_FOLD, 
		(double complex*)dot.m, (double complex*)dot.c);
	return dot;
}

double complex rzdotc(
	int N,
	double complex* x, int incx,
	double complex* y, int incy
) {
	I_double_Complex dot;
	zisetzero(DEFAULT_FOLD, &dot);
	zdotcI1(N, x, incx, y, incy, DEFAULT_FOLD, 
		(double complex*)dot.m, (double complex*)dot.c);
    double complex ret;
    zziconv_sub(DEFAULT_FOLD, &dot, &ret);
	return ret;
}

I_double_Complex zdotuI(
	int N,
	double complex* x, int incx,
	double complex* y, int incy
) {
	I_double_Complex dot;
	zisetzero(DEFAULT_FOLD, &dot);
	zdotuI1(N, x, incx, y, incy, DEFAULT_FOLD, 
		(double complex*)dot.m, (double complex*)dot.c);
	return dot;
}

double complex rzdotu(
	int N,
	double complex* x, int incx,
	double complex* y, int incy
) {
	I_double_Complex dot;
	zisetzero(DEFAULT_FOLD, &dot);
	zdotuI1(N, x, incx, y, incy, DEFAULT_FOLD, 
		(double complex*)dot.m, (double complex*)dot.c);
    double complex ret;
    zziconv_sub(DEFAULT_FOLD, &dot, &ret);
	return ret;
}

