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

void sanity_check(int n, float* x, float* y, int* status) {
	// GENERATE DATA
	int i;
	float small = 1.0 / (1024.0 * 1024.0);			// 2^-20
	float big   = 1024.0 * 1024.0 * 1024.0 * 32;	// 2^35
	scomplex ref[5];
	float res;
	for (i = 0;  i < 5; i++) status[i] = 0;

	// big at the begining
	for (i = 0; i < n; i++) {
		x[4*i] = sin(M_PI * ((float)i / (float)n));
		x[4*i] *= rsign();
		x[4*i + 1] = cos(M_PI * ((float)i / (float)n));
		x[4*i + 1] *= rsign();
	}

	ccheck_reproducibility(n, (scomplex*)x, 2, (scomplex*)y, 1, status,
		(scomplex*)ref);
}

