#include "indexedBLAS.h"
#include <stdlib.h>

void sIAccInit(sIAccum *acc, int BUFFER_SIZE) {
	acc->size = BUFFER_SIZE;
	acc->counter = 0;
	sISetZero(acc->v);

	if (BUFFER_SIZE < 1) BUFFER_SIZE = 128;
	acc->BUFFER = (float*) malloc(BUFFER_SIZE * sizeof(float));
}

void sIAccReset(sIAccum *acc) {
	acc->counter = 0;
	sISetZero(acc->v);
}

void sIAccDestroy(sIAccum *acc) {
	if (acc->BUFFER != NULL) {
		free (acc->BUFFER);
		acc->BUFFER = NULL;
	}
}

void sIAccumulate(sIAccum *acc, float x) {
	int c = acc->counter;
	acc->BUFFER[c] = x;
	c++;
	if (c == acc->size) {
		ssumI1(c, acc->BUFFER, 1, 3, sIWidth(), acc->v.m, acc->v.c);
		c = 0;
	}
	acc->counter = c;
}

void sIAccumulates(sIAccum *acc, int n, float* v, int inc) {
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
		ssumI1(lN + acc->counter, acc->BUFFER, 1,
			0, sIWidth(), acc->v.m, acc->v.c);

		// reset
		acc->counter = 0;
		n -= lN;
		v += lN * inc;
	}
}

float sIAccExtract(sIAccum *acc) {
	int c = acc->counter;
	if (c > 0) {
		ssumI1(c, acc->BUFFER, 1, 3, sIWidth(), acc->v.m, acc->v.c);
	}
	return Iconv2f((acc->v));
}

