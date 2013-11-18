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


void sanity_check_s(int n, float* x, float* y, int* status) {
	// GENERATE DATA
	int i;
	float ref[4];
	for (i = 0;  i < 4; i++) status[i] = 0;

	// big at the begining
	for (i = 0; i < n; i++) {
		x[i] = sin(M_PI * ((double)i / (double)n));
	}

	scheck_reproducibility(n, x, y, status, ref);

}

