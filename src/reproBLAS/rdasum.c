/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

I_double dasumI(int N, double* v, int inc) {
	I_double sum;
	dISetZero(sum);
	dasumI1(N, v, inc, DEFAULT_FOLD, 0, sum.m, sum.c);
	return sum;
}

double rdasum(int N, double* v, int inc) {
	I_double sum;
	dISetZero(sum);
	dasumI1(N, v, inc, DEFAULT_FOLD, 0, sum.m, sum.c);
	return Iconv2d(sum);
}

