#include "idxdBLAS.h"
#include <complex.h>
#include <stdlib.h>

void dIAccInit(dIAccum *acc, int BUFFER_SIZE) {
	acc->size = BUFFER_SIZE;
	acc->counter = 0;
	dISetZero(acc->v);

	if (BUFFER_SIZE < 1) BUFFER_SIZE = 128;
	acc->BUFFER = (double*) malloc(BUFFER_SIZE * sizeof(double));
}

void dIAccReset(dIAccum *acc) {
	acc->counter = 0;
	dISetZero(acc->v);
}

void dIAccDestroy(dIAccum *acc) {
	if (acc->BUFFER != NULL) {
		free (acc->BUFFER);
		acc->BUFFER = NULL;
	}
}

void dIAccumulate(dIAccum *acc, double x) {
	int c = acc->counter;
	acc->BUFFER[c] = x;
	c++;
	if (c == acc->size) {
		dsumI1(c, acc->BUFFER, 1, 3, acc->v.m, acc->v.c);
		c = 0;
	}
	acc->counter = c;
}

void dIAccumulates(dIAccum *acc, int n, double* v, int inc) {
	int lN;
	int i;
	// NAIVE IMPLEMENTATION
	while (n > 0) {
		lN = acc->size - acc->counter;
		lN = lN < n ? lN : n;
		// copy
		for (i = 0; i < lN; i++) {
			acc->BUFFER[i] = v[i * inc];
		}
		dsumI1(lN + acc->counter, acc->BUFFER, 1,
			0, acc->v.m, acc->v.c);

		// reset
		acc->counter = 0;
		n -= lN;
		v += lN * inc;
	}
}

double dIAccExtract(dIAccum *acc) {
	int c = acc->counter;
	if (c > 0) {
		dsumI1(c, acc->BUFFER, 1, 3,  acc->v.m, acc->v.c);
	}
	return idxd_ddiconv(&(acc->v), DEFAULT_FOLD);
}

