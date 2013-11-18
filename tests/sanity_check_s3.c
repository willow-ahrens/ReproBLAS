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

void sanity_check_s(int n, float* x, float* y, int* status) {
	// GENERATE DATA
	int i;
	float small = 1.0 / (1024.0 * 1024.0);			// 2^-20
	float big   = 1024.0 * 1024.0 * 1024.0 * 32;	// 2^35
	float ref[4];
	float res;
	for (i = 0;  i < 4; i++) status[i] = 0;
	float normal = 0.0;
	float anormal = 0.0;

	// big at the begining
	for (i = 0; i < n; i++) {
		x[i] = small + (big - small) * i / n;
//		x[i] *= rsign();

		normal += x[i];
		anormal += fabs(x[i]);
	}
	ref[0] = rssum(n, x, 1);
	ref[1] = rsasum(n, x, 1);
	ref[2] = rsdot(n, x, 1, y, 1);

	// reverse
	sreverse(n, x);

	res = rssum(n, x, 1);
	if (res != ref[0]) status[0]++;

	res = rsasum(n, x, 1);
	if (res != ref[1]) {
		status[1]++;
		printf("\n reverse: %g %g %g %g ", res, ref[1], anormal, res - ref[1]);
	}

	res = rsdot(n, x, 1, y, 1);
	if (res != ref[2]) {
		status[2]++;
		printf("\n reverse: %g %g %g %g ", res, ref[2], normal, res - ref[2]);
	}

	// shuffle
	sshuffle(n, x);

	res = rssum(n, x, 1);
	if (res != ref[0]) status[0]++;

	res = rsasum(n, x, 1);
	if (res != ref[1]) {
		status[1]++;
		printf("\n shuffle: %g %g %g %g ", res, ref[1], anormal, res - ref[1]);
	}

	res = rsdot(n, x, 1, y, 1);
	if (res != ref[2]) {
		status[2]++;
		printf("\n shuffle %g %g %g %g ", res, ref[2], normal, res - ref[2]);
	}
}

