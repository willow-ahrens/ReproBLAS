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

void sanity_check(int n, double* x, double* y, int* status) {
	// GENERATE DATA
	int i;
	double small = 1.0 / (1024.0 * 1024.0);			// 2^-20
	double big   = 1024.0 * 1024.0 * 1024.0 * 32;	// 2^35
	double ref   = (n - 2) * small;
	double refa  = ((n - 2) * small) + 2 * big;
	double res;
	for (i = 0; i < n; i++) {
		x[i] = small;
	}
	for (i = 0;  i < 4; i++) status[i] = 0;

	// big at the begining
	x[0]   = big;
	x[n-1] = -big;

	res = rdsum(n, x, 1);
	if (res != ref) status[0]++;

	res = rdasum(n, x, 1);
	if (res != refa) status[1]++;

	res = rddot(n, x, 1, y, 1);
	if (res != ref) status[2]++;

	x[0] = small;
	x[n-1] = small;

	// big at the end
	x[n/2] = big;
	x[n-1] = -big;

	res = rdsum(n, x, 1);
	if (res != ref) status[0]++;

	res = rdasum(n, x, 1);
	if (res != refa) status[1]++;

	res = rddot(n, x, 1, y, 1);
	if (res != ref) status[2]++;

	x[n/2] = small;
	x[n-1] = small;

	// big in the middle
	x[0] = big;
	x[n/2] = -big;

	res = rdsum(n, x, 1);
	if (res != ref) status[0]++;

	res = rdasum(n, x, 1);
	if (res != refa) status[1]++;

	res = rddot(n, x, 1, y, 1);
	if (res != ref) status[2]++;

	x[0] = small;
	x[n/2] = small;

}

