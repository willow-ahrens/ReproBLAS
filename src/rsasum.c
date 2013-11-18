/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rblas1.h"

I_float sasumI(int N, float* v, int inc) {
	I_float sum;

	sISetZero(sum);
	sasumI1(N, v, inc, DEFAULT_FOLD, 0, sum.m, sum.c);

	return sum;
}

float rsasum(int N, float* v, int inc) {
	I_float sum;

	sISetZero(sum);
	sasumI1(N, v, inc, DEFAULT_FOLD, 0, sum.m, sum.c);

	return Iconv2f(sum);
}

