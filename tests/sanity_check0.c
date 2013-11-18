#include <rblas.h>
#include <stdio.h>
#include <stdlib.h>

char* tname = "TEST 0";
char* tdesc = "Unit in the first place";

const char* name() {
	return tname;
}
const char* desc() {
	return tdesc;
}

void sanity_check(int n, double* x, double* y, int* status) {
	// GENERATE DATA
	int i;
	n = 10;
	for (i = 0; i < 4; i++) status[i] = 0;
	double small = 1.0 / (1024.0 * 1024.0);			// 2^-20
	double bin = 1.0;
	double ufpx;
	for (i = 0; i < n; i++) {
		x[i] = 3 * bin;
		bin *= 2;

		// check
		ufpx = ufp(x[i]);
		if (ufpx != bin) {
			status[0] = 1;
			printf("UFP(%g) = %g != %g\n", x[i], ufpx, bin);
		}
	}
}

