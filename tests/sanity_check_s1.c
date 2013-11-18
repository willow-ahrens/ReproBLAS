#include <rblas.h>
#include <stdio.h>
#include <stdlib.h>

char* tname = "TEST 1";
char* tdesc = "1 Big number";

const char* name() {
	return tname;
}
const char* desc() {
	return tdesc;
}

void sanity_check_s(int n, float* x, float* y, int* status) {
	// GENERATE DATA
	int i;
	float small = 1.0 / 1024.0;			// 2^-10
	float big   = 1024.0 * 32;		// 2^15
	float ref   = (n - 1) * small;
	float res;
	ref += big;
	for (i = 0; i < n; i++) {
		x[i] = small;
	}
	for (i = 0;  i < 4; i++) status[i] = 0;

	// big at the begining
	x[0] = big;

	res = rssum(n, x, 1);
	if (res != ref) {
		printf("%g <> %g \n", ref, res);
		status[0]++;
	}
	
	res = rsasum(n, x, 1);
	if (res != ref) {
		printf("%g <> %g \n", ref, res);
		status[1]++;
	}

	res = rsdot(n, x, 1, y, 1);
	if (res != ref) status[2]++;

	x[0] = small;

	// big at the end
	x[n-1] = big;

	res = rssum(n, x, 1);
	if (res != ref) status[0]++;

	res = rsasum(n, x, 1);
	if (res != ref) status[1]++;

	res = rsdot(n, x, 1, y, 1);
	if (res != ref) status[2]++;

	x[n-1] = small;

	// big in the middle
	x[n/2] = big;

	res = rssum(n, x, 1);
	if (res != ref) status[0]++;

	res = rsasum(n, x, 1);
	if (res != ref) status[1]++;

	res = rsdot(n, x, 1, y, 1);
	if (res != ref) status[2]++;

	x[n/2] = small;

}

