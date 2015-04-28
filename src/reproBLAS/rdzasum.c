/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "reproBLAS.h"
#include "indexedBLAS.h"

I_double dzasumI(int N, double complex* v, int inc) {
	I_double sum;
	disetzero(DEFAULT_FOLD, &sum);
	dzasumI1(N, v, inc, DEFAULT_FOLD, sum.m, sum.c);
	return sum;
}

double rdzasum(int N, double complex* v, int inc) {
	I_double sum;
	disetzero(DEFAULT_FOLD, &sum);
	dzasumI1(N, v, inc, DEFAULT_FOLD,  sum.m, sum.c);
	return ddiconv(&sum, DEFAULT_FOLD);
}

