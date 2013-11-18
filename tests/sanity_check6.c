#include <rblas.h>
#include <stdio.h>
#include <stdlib.h>
#include "debug.h"
#include <math.h>

char* tname = "TEST 6";
//char* tdesc = "Condition Number = 1e8, 1e12, 1e15";
char tdesc[256];

const char* name() {
	return tname;
}
const char* desc() {
	return tdesc;
}


void sanity_check(int n, double* x, double* y, int* status) {
	// GENERATE DATA
	int i;
	double ref[4];
	double res;
	for (i = 0;  i < 4; i++) status[i] = 0;

	drandomseed();

	sprintf(tdesc, "RCOND :");

	x[0] = 1e8;
	dgenvec(n, x, 5, 1.0);

	dcheck_reproducibility(n, x, y, status, ref);
	sprintf(tdesc, "%s, %g", tdesc, fabs(ref[0]) /  ref[1]);

	x[0] = 1e12;
	dgenvec(n, x, 5, 1.0);

	dcheck_reproducibility(n, x, y, status, ref);
	sprintf(tdesc, "%s, %g", tdesc, fabs(ref[0]) /  ref[1]);

	x[0] = 1e15;
	dgenvec(n, x, 5, 1.0);

	dcheck_reproducibility(n, x, y, status, ref);
	sprintf(tdesc, "%s, %g", tdesc, fabs(ref[0]) /  ref[1]);
}

