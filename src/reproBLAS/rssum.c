/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "IndexedBLAS.h"

I_float ssumI(int N, float* v, int inc) {
	I_float sum;

	sISetZero(sum);
	ssumI1(N, v, inc, DEFAULT_FOLD, 0, sum.m, sum.c);

	return sum;
}

float rssum(int N, float* v, int inc) {
	I_float sum;

	sISetZero(sum);
	ssumI1(N, v, inc, DEFAULT_FOLD, 0, sum.m, sum.c);

	return Iconv2f(sum);
}

