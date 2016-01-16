#include "idxdBLAS.h"
#include <stdlib.h>

void zIAccInit(zIAccum *acc, int BUFFER_SIZE) {
	acc->size = BUFFER_SIZE;
	acc->counter = 0;
	zISetZero(acc->v);

	if (BUFFER_SIZE < 1) BUFFER_SIZE = 128;
	acc->BUFFER = (double complex*) malloc(BUFFER_SIZE * sizeof(double complex));
}

void zIAccReset(zIAccum *acc) {
	acc->counter = 0;
	zISetZero(acc->v);
}

void zIAccDestroy(zIAccum *acc) {
	if (acc->BUFFER != NULL) {
		free (acc->BUFFER);
		acc->BUFFER = NULL;
	}
}

void zIAccumulate(zIAccum *acc, double complex x) {
	int c = acc->counter;
	acc->BUFFER[c] = x;
	c++;
	if (c == acc->size) {
		zsumI1(c, acc->BUFFER, 1, 0, 
			(double complex*)acc->v.m, (double complex*)acc->v.c);
		c = 0;
	}
	acc->counter = c;
}

void zIAccumulates(zIAccum *acc, int n, double complex* v, int inc) {
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
		zsumI1(lN + acc->counter, acc->BUFFER, 1,
			0, 
			(double complex*)acc->v.m, (double complex*)acc->v.c);

		// reset
		acc->counter = 0;
		n -= lN;
		v += lN * inc;
	}
}

double complex zIAccExtract(zIAccum *acc) {
	int c = acc->counter;
	if (c > 0) {
		zsumI1(c, acc->BUFFER, 1, 0, 
			(double complex*)acc->v.m, (double complex*)acc->v.c);
	}
	return Iconv2z((acc->v));
}

