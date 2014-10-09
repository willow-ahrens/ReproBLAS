/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rblas1.h"

I_float_Complex csumI(int N, float complex* v, int inc) {
	I_float_Complex sum;
	cISetZero(sum);
	csumI1(N, v, inc, DEFAULT_FOLD, 0, (float complex*)sum.m, sum.c);
	return sum;
}

float complex rcsum(int N, float complex* v, int inc) {
	I_float_Complex sum;
	cISetZero(sum);
	csumI1(N, v, inc, DEFAULT_FOLD, 0, (float complex*)sum.m, sum.c);
	return Iconv2c(sum);
}

