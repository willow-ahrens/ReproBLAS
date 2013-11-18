#include <rblas.h>
#include <stdio.h>
#include <stdlib.h>
#include "debug.h"
#include <math.h>

char* tname = "TEST 6";
char* tdesc = "cond(real), cond(imaginary) = 1e8, 1e12, 1e15";

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
	double res;
	for (i = 0;  i < 5; i++) status[i] = 0;

	drandomseed();

	y[0] = 1e8;
	y[n] = 1e8;
	dgenvec(n, y, 5, 1.0);
	dgenvec(n, y + n, 5, 1.0);
	for (i = 0; i < n; i++) {
		x[4*i] = y[i];
		x[4*i + 1] = y[n+i];
	}
	for (i = 0; i < 2 *n; i++) y[i] = 1.0;

	zcheck_reproducibility(n, (dcomplex*)x, 2, (dcomplex*)y, 1, status,
		(dcomplex*)ref);

	y[0] = 1e12;
	y[n] = 1e12;
	dgenvec(n, y, 5, 1.0);
	dgenvec(n, y + n, 5, 1.0);
	// big at the begining
	for (i = 0; i < n; i++) {
		x[4*i] = y[i];
		x[4*i + 1] = y[n+i];
	}
	for (i = 0; i < 2 *n; i++) y[i] = 1.0;

	zcheck_reproducibility(n, (dcomplex*)x, 2, (dcomplex*)y, 1, status,
		(dcomplex*)ref);

	y[0] = 1e15;
	y[n] = 1e15;
	dgenvec(n, y, 5, 1.0);
	dgenvec(n, y + n, 5, 1.0);
	// big at the begining
	for (i = 0; i < n; i++) {
		x[4*i] = y[i];
		x[4*i + 1] = y[n+i];
	}
	for (i = 0; i < 2 *n; i++) y[i] = 1.0;

	zcheck_reproducibility(n, (dcomplex*)x, 2, (dcomplex*)y, 1, status,
		(dcomplex*)ref);
}

