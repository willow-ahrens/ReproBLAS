#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rblas1.h"

I_double_Complex zsumI(int N, double complex* v, int inc) {
	I_double_Complex sum;
	zISetZero(sum);
	zsumI1(N, v, inc, DEFAULT_FOLD, 0, (double complex*)sum.m, (double complex*)sum.c);
	return sum;
}

double complex rzsum(int N, double complex* v, int inc) {
	I_double_Complex sum;
	zISetZero(sum);
	zsumI1(N, v, inc, DEFAULT_FOLD, 0, (double complex*)sum.m, (double complex*)sum.c);
	return Iconv2z(sum);
}

