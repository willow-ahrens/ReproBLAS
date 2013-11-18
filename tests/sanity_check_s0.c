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

void sanity_check_s(int n, float* x, float* y, int* status) {
	// GENERATE DATA
	int i;
	n = 10;
	for (i = 0; i < 4; i++) status[i] = 0;
	float small = 1.0 / (1024.0 * 1024.0);			// 2^-20
	float bin = 1.0;
	float ufp;
//	printf("size of I_float: %ld \n", sizeof(I_float));
	for (i = 0; i < n; i++) {
		x[i] = 3 * bin;
		bin *= 2;

		// check
		ufp = ufpf(x[i]);
		if (ufp != bin) {
			status[0] = 1;
			printf("UFP(%g) = %g != %g\n", x[i], ufp, bin);
		}
	}
}

