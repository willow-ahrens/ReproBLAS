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


void sanity_check_s(int n, float* x, float* y, int* status) {
	// GENERATE DATA
	int i;
	float ref[4];
	float res;
	for (i = 0;  i < 4; i++) status[i] = 0;

	sprintf(tdesc, "RCOND :");

	x[0] = 1e3;
	sgenvec(n, x, 5, 1.0);

	scheck_reproducibility(n, x, y, status, ref);
	sprintf(tdesc, "%s, %g", tdesc, fabs(ref[0]) /  ref[1]);

	x[0] = 1e5;
	sgenvec(n, x, 5, 1.0);

	scheck_reproducibility(n, x, y, status, ref);
	sprintf(tdesc, "%s, %g", tdesc, fabs(ref[0]) /  ref[1]);

	x[0] = 1e8;
	sgenvec(n, x, 5, 1.0);

	scheck_reproducibility(n, x, y, status, ref);
	sprintf(tdesc, "%s, %g", tdesc, fabs(ref[0]) /  ref[1]);
}

