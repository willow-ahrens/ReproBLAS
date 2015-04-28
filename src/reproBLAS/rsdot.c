#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

I_float sdotI(int N, float* x, int incx, float* y, int incy) {
	I_float dot;

	sisetzero(DEFAULT_FOLD, &dot);
	sdotI1(N, x, incx, y, incy, DEFAULT_FOLD, dot.m, dot.c);

	return dot;
}

float rsdot(int N, float* x, int incx, float* y, int incy) {
	I_float dot;

	sisetzero(DEFAULT_FOLD, &dot);
	sdotI1(N, x, incx, y, incy, DEFAULT_FOLD, dot.m, dot.c);

	return ssiconv(DEFAULT_FOLD, &dot);
}

