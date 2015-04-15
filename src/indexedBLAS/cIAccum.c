#include "indexedBLAS.h"
#include <stdlib.h>

void cIAccInit(cIAccum *acc, int BUFFER_SIZE) {
	acc->size = BUFFER_SIZE;
	acc->counter = 0;
	cISetZero(acc->v);

	if (BUFFER_SIZE < 1) BUFFER_SIZE = 128;
	acc->BUFFER = (float complex*) malloc(BUFFER_SIZE * sizeof(float complex));
}

void cIAccReset(cIAccum *acc) {
	acc->counter = 0;
	cISetZero(acc->v);
}

void cIAccDestroy(cIAccum *acc) {
	if (acc->BUFFER != NULL) {
		free (acc->BUFFER);
		acc->BUFFER = NULL;
	}
}

void cIAccumulate(cIAccum *acc, float complex x) {
	int c = acc->counter;
	acc->BUFFER[c] = x;
	c++;
	if (c == acc->size) {
		csumI1(c, acc->BUFFER, 1, 0, 
			(float complex*)acc->v.m, acc->v.c);
		c = 0;
	}
	acc->counter = c;
}

void cIAccumulates(cIAccum *acc, int n, float complex* v, int inc) {
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
		csumI1(lN + acc->counter, acc->BUFFER, 1,
			0, 
			(float complex*)acc->v.m, acc->v.c);

		// reset
		acc->counter = 0;
		n -= lN;
		v += lN * inc;
	}
}

float complex cIAccExtract(cIAccum *acc) {
	int c = acc->counter;
	if (c > 0) {
		csumI1(c, acc->BUFFER, 1, 0, 
			(float complex*)acc->v.m, acc->v.c);
	}
	return Iconv2c((acc->v));
}

