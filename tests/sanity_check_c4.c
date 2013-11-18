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

void sanity_check(int n, float* x, float* y, int* status) {
	// GENERATE DATA
	int i;
	float small = 1.0 / 1024.0 ;			// 2^-10
	float big   = 1024.0 * 1024.0;	// 2^20
	dcomplex ref[5];
	float res;
	for (i = 0;  i < 5; i++) status[i] = 0;

	// big at the begining
	int inc = 2;
	for (i = 0; i < n; i++) {
		x[2*i*inc] = small + (big - small) * drand48();
		x[2*i*inc] *= rsign();
		x[2*i*inc + 1] = small + (big - small) * drand48();
		x[2*i*inc + 1] *= rsign();
//		x[4*i + 2] = x[4*i+3] = 0.0;
	}

	ccheck_reproducibility(n, (scomplex*)x, inc, (scomplex*)y, 1, status, (scomplex*)ref);
}

