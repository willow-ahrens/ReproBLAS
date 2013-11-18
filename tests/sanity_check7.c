#include <rblas.h>
#include <stdio.h>
#include <stdlib.h>
#include "debug.h"
#include <math.h>

char* tname = "TEST 7";
char* tdesc = "x[i] = sin(PI * i / n)";

const char* name() {
	return tname;
}
const char* desc() {
	return tdesc;
}


void sanity_check(int n, double* x, double* y, int* status) {
	// GENERATE DATA
	int i;
	double small = 1.0 / (1024.0 * 1024.0);			// 2^-20
	double big   = 1024.0 * 1024.0 * 1024.0 * 32;	// 2^35
	double ref[4];
	double res;
	for (i = 0;  i < 4; i++) status[i] = 0;

	// big at the begining
	double t1, t2, x1;
	for (i = 0; i < n; i++) {
		x[i] = sin(M_PI * ((double)i / (double)n));
	}

	dcheck_reproducibility(n, x, y, status, ref);

}

