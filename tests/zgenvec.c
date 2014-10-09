#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "debug.h"
#include "../src/rblas1.h"

#define SWAP_(A,B,T) { T = A; A = B; B = T; }

#ifdef _MATH__COMPLEX_
#define fabsz(X) cabs(X)
#else
#define fabsz(X) (fabs(X.real)+fabs(X.imag))
#endif

void zreverse(int N, dcomplex* v, int inc) {
	int i;
	dcomplex t;
	for (i = 0 ; i < N / 2; i++)
	{
		SWAP_(v[i*inc], v[(N -i-1) * inc], t);	
	}
}

void zshuffle(int N, dcomplex* v, int inc) {
	int i;
	dcomplex t;
	int gap;
	for (i = 0 ; i < N - 2; i++)
	{
		gap = rand() % (N - i);
		SWAP_(v[i * inc], v[(i + gap) * inc], t);	
	}
}

// SLOW SORTING
//    ord: 1  INCREASING
//         -1 DECREASING
//         2  INCREASING MAGNITUDE
//         -2 DECREASING MAGNITUDE
void zsort_bubble(int N, dcomplex* v, int inc, int ord) {
	int conted = 1;
	int i;
	dcomplex t;

	while (conted == 1) {
		conted = 0;

		switch (ord) {
		case 1:
			for (i = 0; i < N-1; i++) {
				if (ZREAL_(v[i*inc]) > ZREAL_(v[i*inc+inc])) {
					conted = 1;
					SWAP_(v[i*inc], v[i*inc+inc], t);
				}
			}
			break;
		case -1:
			for (i = 0; i < N-1; i++) {
				if (ZREAL_(v[i*inc]) < ZREAL_(v[i*inc+inc])) {
					conted = 1;
					SWAP_(v[i*inc], v[i*inc+inc], t);
				}
			}
			break;
		case 2:
			for (i = 0; i < N-1; i++) {
				if (fabsz(v[i*inc]) > fabsz(v[i*inc+inc])) {
					conted = 1;
					SWAP_(v[i*inc], v[i*inc+inc], t);
				}
			}
			break;
		case -2:
			for (i = 0; i < N-1; i++) {
				if (fabsz(v[i*inc]) < fabsz(v[i*inc+inc])) {
					conted = 1;
					SWAP_(v[i*inc], v[i*inc+inc], t);
				}
			}
			break;
		}
	}
}

// merge v1, v2 into v
void zmerge_(
	int N,
	int N1, dcomplex* v1, int inc1, int* ptr_id1,
	int N2, dcomplex* v2, int inc2, int* ptr_id2,
	dcomplex* v, int inc,
	int ord
) {
	int id1, id2, id;

	id1 = ptr_id1[0];
	id2 = ptr_id2[0];

	switch (ord) {
		case 1:
			for (id = 0; id < N && id1 < N1 && id2 < N2; id++) {
				if (ZREAL_(v1[id1*inc1]) < ZREAL_(v2[id2*inc2])) {
					v[id*inc] = v1[inc1*id1++];
				}
				else
					v[id*inc] = v2[inc2*id2++];
			}
			break;
		case -1:
			for (id = 0; id < N && id1 < N1 && id2 < N2; id++) {
				if (ZREAL_(v1[id1*inc]) > ZREAL_(v2[id2*inc])) {
					v[id*inc] = v1[inc1*id1++];
				}
				else
					v1[id] = v2[id2++];
			}
			break;
		case 2:
			for (id = 0; id < N && id1 < N1 && id2 < N2; id++) {
				if (fabsz(v1[id1*inc1]) < fabsz(v2[id2*inc2])) {
					v[id*inc] = v1[inc1*id1++];
				}
				else
					v[id*inc] = v2[inc2*id2++];
			}
			break;
		case -2:
			for (id = 0; id < N && id1 < N1 && id2 < N2; id++) {
				if (fabsz(v1[id1*inc1]) > fabsz(v2[id2*inc2])) {
					v[id*inc] = v1[inc1*id1++];
				}
				else
					v[id*inc] = v2[inc2*id2++];
			}
			break;
	}

	if (id1 == N1 && id < N) {
		if (v != v2) {
			for (; id < N; id++, id2++)
				v[id*inc] = v2[id2 * inc2];
		}
		id2 += N - id;
	}
	else if (id2 == N2 && id < N) {
		if (v != v1) {
			for (; id < N; id++, id1++)
				v[id*inc] = v1[id1 * inc1];
		}
		id1 += N - id;
	}

	ptr_id1[0] = id1;
	ptr_id2[0] = id2;
}

void zmerge_sort_(
	int N1, dcomplex* v1, int inc1,
	int N2, dcomplex* v2, int inc2,
	int ord,
	dcomplex* buffer
) {
	int id1, id2;
	int i;

	id1 = id2 = 0;
	//if (N1 <= N2) {
		// COPY v1 to buffer
		for (i = 0; i < N1; i++)
			buffer[i] = v1[i * inc1];

		// FIRST PHASE: MERGE BUFFER,V2 TO V1
		zmerge_(N1, N1, buffer, 1, &id1, N2, v2, inc2, &id2, v1, inc1, ord);

		// buffer empty: done
		if (id1 == N1)
			return;

		// V2 empty: copy buffer to v2
		if (id2 == N2) {
			for (i = 0; i < N2; i++)
				v2[i * inc2] = buffer[i];
			return;
		}

		//SECOND PHASE: MERGE BUFFER+ID1, V2+ID2 TO V2
		zmerge_(N2, N1, buffer, 1, &id1, N2, v2, inc2, &id2, v2, inc2, ord);

		return;
	//}
}

void zsort_merge_(
	int N,
	dcomplex* v, int inc,
	int ord,
	dcomplex* buffer
) {
	if (N <= 32) {
		zsort_bubble(N, v, inc, ord);
		return;
	}

	int mid = N /2;

	// left recursive call
	zsort_merge_(mid, v, inc, ord, buffer);
	// right recursive call
	zsort_merge_(N-mid, v+mid*inc, inc, ord, buffer);
	// merge
	zmerge_sort_(mid, v, inc, N-mid, v+mid*inc, inc, ord, buffer);
}

void zsort_merge(
	int N,
	dcomplex* v, int inc,
	int ord
) {
	dcomplex* buffer = (dcomplex*) malloc((N/2) * sizeof(dcomplex));

	zsort_merge_(N, v, inc, ord, buffer);

	free(buffer);
}

