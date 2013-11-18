#include <rblas.h>
#include <stdio.h>
#include <stdlib.h>
#include "debug.h"

char* tname = "TEST 4";
char* tdesc = "small + (big-small) * drand48()";

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
	dcomplex ref[5];
	double res;
	for (i = 0;  i < 5; i++) status[i] = 0;

	// big at the begining
	for (i = 0; i < n; i++) {
		x[4*i] = small + (big - small) * drand48();
		x[4*i] *= rsign();
		x[4*i + 1] = small + (big - small) * drand48();
		x[4*i + 1] *= rsign();
	}

	zcheck_reproducibility(n, (dcomplex*)x, 2, (dcomplex*)y, 1, status,
		(dcomplex*)ref);
}

