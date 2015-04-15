/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

#ifdef __SSE__
#	include <xmmintrin.h>
#endif

void sIRenorm1(int n, float* X, float* leading, int inc) {
	int i;
	float M;
	float x;
	for (i = 0; i < n; i++, X+=inc, leading+=inc) {
		x = X[0];
		if (x == 0.0)
			continue;

		M = ufpf(x);
		if (x >= (M * 1.75)) {
			X[0] -= M * 0.25;
			leading[0] += 1;
		}
		else if (x < (M * 1.25)) {
			X[0] += M * 0.5;
			leading[0] -= 2;
		}
		else if (x < (M * 1.5)) {
			X[0] += M * 0.25;
			leading[0] -= 1;
		}
	}
}

void cIRenorm1(int K, float complex* X, float* C, int INC) {
	sIRenorm1(K, (float*)X,     C,     2*INC);
	sIRenorm1(K, (float*)X + 1, C + 1, 2*INC);
}

