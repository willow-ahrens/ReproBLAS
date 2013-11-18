#include <rblas.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "debug.h"

char* tname = "TEST 7";
char* tdesc = "sin(PI * i /n) + I * cos(PI * i / n)";

const char* name() {
	return tname;
}
const char* desc() {
	return tdesc;
}

void sanity_check(int n, double* x, double* y, int* status) {
	// GENERATE DATA
	int i;
	dcomplex ref[5];
	for (i = 0;  i < 5; i++) status[i] = 0;

	// big at the begining
	for (i = 0; i < n; i++) {
		x[4*i] = sin(M_PI * ((double)i / (double)n));
		x[4*i] *= rsign();
		x[4*i + 1] = cos(M_PI * ((double)i / (double)n));
		x[4*i + 1] *= rsign();
	}

	zcheck_reproducibility(n, (dcomplex*)x, 2, (dcomplex*)y, 1, status,
		(dcomplex*)ref);
}

