#include <rblas.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "debug.h"

char* tname = "TEST 3";
char* tdesc = "small + (big-small)*i/n";

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
	double ref[4];
	double res;
	for (i = 0;  i < 4; i++) status[i] = 0;
	double normal = 0.0;
	double anormal = 0.0;

	// big at the begining
	for (i = 0; i < n; i++) {
		x[i] = small + (big - small) * i / n;
//		x[i] *= rsign();

		normal += x[i];
		anormal += fabs(x[i]);
	}
	ref[0] = rdsum(n, x, 1);
	ref[1] = rdasum(n, x, 1);
	ref[2] = rddot(n, x, 1, y, 1);

	// reverse
	dreverse(n, x);

	res = rdsum(n, x, 1);
	if (res != ref[0]) status[0]++;

	res = rdasum(n, x, 1);
	if (res != ref[1]) {
		status[1]++;
		printf("\n reverse: %g %g %g %g ", res, ref[1], anormal, res - ref[1]);
	}

	res = rddot(n, x, 1, y, 1);
	if (res != ref[2]) {
		status[2]++;
		printf("\n reverse: %g %g %g %g ", res, ref[2], normal, res - ref[2]);
	}

	// shuffle
	dshuffle(n, x);

	res = rdsum(n, x, 1);
	if (res != ref[0]) status[0]++;

	res = rdasum(n, x, 1);
	if (res != ref[1]) {
		status[1]++;
		printf("\n shuffle: %g %g %g %g ", res, ref[1], anormal, res - ref[1]);
	}

	res = rddot(n, x, 1, y, 1);
	if (res != ref[2]) {
		status[2]++;
		printf("\n shuffle %g %g %g %g ", res, ref[2], normal, res - ref[2]);
	}
}

