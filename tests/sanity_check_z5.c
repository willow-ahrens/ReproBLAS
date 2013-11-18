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


void sanity_check(int n, double* x, double* y, int* status) {
	// GENERATE DATA
	int i;
	dcomplex ref[5];
	for (i = 0;  i < 5; i++) status[i] = 0;

	// big at the begining
	double t1, t2, x1;
	for (i = 0; i < 4*n; i++) {
		t1 = drand48();
		t2 = drand48();
		x1 = sqrt(-2.0 * log(t1)) * cos(2.0 * M_PI * t2);
		x[i] = x1 * rsign();
	}

	zcheck_reproducibility(n, (dcomplex*)x, 2, (dcomplex*)y, 1, status,
		(dcomplex*)ref);

}

