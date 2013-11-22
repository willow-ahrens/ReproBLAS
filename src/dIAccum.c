#include "rblas1.h"
#include "Iblas1.h"
#include <stdlib.h>

void dInitAccum(dIAccum *acc, int BUFFER_SIZE) {
	acc->size = BUFFER_SIZE;
	acc->counter = 0;
	dISetZero(acc->v);

	if (BUFFER_SIZE < 1) BUFFER_SIZE = 128;
	acc->BUFFER = (double*) malloc(BUFFER_SIZE * sizeof(double));
}

void dResetAcc(dIAccum *acc) {
	acc->counter = 0;
	dISetZero(acc->v);
}

void dDestroyAcc(dIAccum *acc) {
	if (acc->BUFFER != NULL) {
		free (acc->BUFFER);
		acc->BUFFER = NULL;
	}
}

void dAccumulate(dIAccum *acc, double x) {
	int c = acc->counter;
	acc->BUFFER[c] = x;
	c++;
	if (c == acc->size) {
		dsumI1(c, acc->BUFFER, 1, 3, dIWidth(), acc->v.m, acc->v.c);
		c = 0;
	}
	acc->counter = c;
}

void dAccumulates(dIAccum *acc, int n, double* v, int inc) {
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
			0, dIWidth(), acc->v.m, acc->v.c);

		// reset
		acc->counter = 0;
		n -= lN;
		v += lN * inc;
	}
}

double dExtractAcc(dIAccum *acc) {
	int c = acc->counter;
	if (c > 0) {
		dsumI1(c, acc->BUFFER, 1, 3, dIWidth(), acc->v.m, acc->v.c);
	}
	return Iconv2d((acc->v));
}

