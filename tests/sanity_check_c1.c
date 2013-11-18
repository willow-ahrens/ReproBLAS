#include <rblas.h>
#include <stdio.h>
#include <stdlib.h>
#include "debug.h"

char* tname = "TEST 1";
char* tdesc = "x[i] = 1";

const char* name() {
	return tname;
}
const char* desc() {
	return tdesc;
}

void sanity_check(int n, float* x, float* y, int* status) {
	// GENERATE DATA
	int i;
	scomplex ref[5];
	float res;
	for (i = 0;  i < 5; i++) status[i] = 0;

	// big at the begining
	int inc = 2;
	for (i = 0; i < n; i++) {
		x[2*inc*i] = 1;
		x[2*inc*i + 1] = x[2*inc*i];
	}

	ccheck_reproducibility(n, (scomplex*)x, inc, (scomplex*)y, 1, status, (scomplex*)ref);
}

