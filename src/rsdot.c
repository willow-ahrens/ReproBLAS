#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rblas1.h"

I_float sdotI(int N, float* x, int incx, float* y, int incy) {
	I_float dot;

	sISetZero(dot);
	sdotI1(N, x, incx, y, incy, DEFAULT_FOLD, 0, dot.m, dot.c);

	return dot;
}

float rsdot(int N, float* x, int incx, float* y, int incy) {
	I_float dot;

	sISetZero(dot);
	sdotI1(N, x, incx, y, incy, DEFAULT_FOLD, 0, dot.m, dot.c);

	return Iconv2f(dot);
}

