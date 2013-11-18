#include <rblas.h>
#include <stdio.h>
#include <stdlib.h>
#include "debug.h"
#include <math.h>

char* tname = "TEST 5";
char* tdesc = "Normal distribution";

const char* name() {
	return tname;
}
const char* desc() {
	return tdesc;
}


void sanity_check_s(int n, float* x, float* y, int* status) {
	// GENERATE DATA
	int i;
	float small = 1.0 / (1024.0 * 1024.0);			// 2^-20
	float big   = 1024.0 * 1024.0 * 1024.0 * 32;	// 2^35
	float ref[4];
	float res;
	for (i = 0;  i < 4; i++) status[i] = 0;

	// big at the begining
	float t1, t2, x1;
	for (i = 0; i < n; i++) {
		t1 = drand48();
		t2 = drand48();
		x1 = sqrt(-2.0 * log(t1)) * cos(2.0 * M_PI * t2);
		x[i] = x1 * rsign();
	}

	scheck_reproducibility(n, x, y, status, ref);

}

