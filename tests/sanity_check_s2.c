#include <rblas.h>
#include <stdio.h>
#include <stdlib.h>

char* tname = "TEST 2";
char* tdesc = "1 Big pos, 1 Big neg";

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
	float ref   = (n - 2) * small;
	float refa  = ((n - 2) * small) + 2 * big;
	float res;
	for (i = 0; i < n; i++) {
		x[i] = small;
	}
	for (i = 0;  i < 4; i++) status[i] = 0;

	// big at the begining
	x[0]   = big;
	x[n-1] = -big;

	ref = rssum(n, x, 1);

	res = rsasum(n, x, 1);
	if (res != refa) status[1]++;

	res = rsdot(n, x, 1, y, 1);
	if (res != ref) status[2]++;

	x[0] = small;
	x[n-1] = small;

	// big at the end
	x[n/2] = big;
	x[n-1] = -big;

	res = rssum(n, x, 1);
	if (res != ref) {
		status[0]++;
		printf("\n n/2,n: %g - %g = %g", res, ref, res - ref);
	}

	res = rsasum(n, x, 1);
	if (res != refa) status[1]++;

	res = rsdot(n, x, 1, y, 1);
	if (res != ref) status[2]++;

	x[n/2] = small;
	x[n-1] = small;

	// big in the middle
	x[0] = big;
	x[n/2] = -big;

	res = rssum(n, x, 1);
	if (res != ref) {
		status[0]++;
		printf("\n 0,n/2: %g - %g = %g", res, ref, res - ref);
	}

	res = rsasum(n, x, 1);
	if (res != refa) status[1]++;

	res = rsdot(n, x, 1, y, 1);
	if (res != ref) status[2]++;

	x[0] = small;
	x[n/2] = small;

}

