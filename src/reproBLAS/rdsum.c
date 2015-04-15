/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

I_double dsumI(int N, double* v, int inc) {
	I_double sum;
	dISetZero(sum);
	dsumI1(N, v, inc, DEFAULT_FOLD, sum.m, sum.c);
	return sum;
}

double rdsum(int N, double* v, int inc) {
	I_double sum;
	dISetZero(sum);
	dsumI1(N, v, inc, DEFAULT_FOLD, sum.m, sum.c);
	return Iconv2d(sum);
}

