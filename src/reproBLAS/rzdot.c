#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "IndexedBLAS.h"

I_double_Complex zdotcI(
	int N,
	double complex* x, int incx,
	double complex* y, int incy
) {
	I_double_Complex dot;
	zISetZero(dot);
	zdotcI1(N, x, incx, y, incy, DEFAULT_FOLD, 0,
		(double complex*)dot.m, (double complex*)dot.c);
	return dot;
}

double complex rzdotc(
	int N,
	double complex* x, int incx,
	double complex* y, int incy
) {
	I_double_Complex dot;
	zISetZero(dot);
	zdotcI1(N, x, incx, y, incy, DEFAULT_FOLD, 0,
		(double complex*)dot.m, (double complex*)dot.c);
	return Iconv2z(dot);
}

I_double_Complex zdotuI(
	int N,
	double complex* x, int incx,
	double complex* y, int incy
) {
	I_double_Complex dot;
	zISetZero(dot);
	zdotuI1(N, x, incx, y, incy, DEFAULT_FOLD, 0,
		(double complex*)dot.m, (double complex*)dot.c);
	return dot;
}

double complex rzdotu(
	int N,
	double complex* x, int incx,
	double complex* y, int incy
) {
	I_double_Complex dot;
	zISetZero(dot);
	zdotuI1(N, x, incx, y, incy, DEFAULT_FOLD, 0,
		(double complex*)dot.m, (double complex*)dot.c);
	return Iconv2z(dot);
}

