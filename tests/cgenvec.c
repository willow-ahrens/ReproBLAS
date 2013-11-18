#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "debug.h"
#include "../src/rblas1.h"

#define SWAP_(A,B,T) { T = A; A = B; B = T; }
#ifdef _MATH__COMPLEX_
#define fabsz(X) cabsf(X)
#else
#define fabsz(X) (fabs(X.real)+fabs(X.imag))
#endif

#define GT_(X,Y) CREAL_(X) > CREAL_(Y)

void creverse(int N, float complex* v, int inc) {
	int i;
	float complex t;
	for (i = 0 ; i < N / 2; i++)
	{
		SWAP_(v[i*inc], v[(N -i-1) * inc], t);	
	}
}

void cshuffle(int N, float complex* v, int inc) {
	int i;
	float complex t;
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
void csort_bubble(int N, float complex* v, int inc, int ord) {
	int conted = 1;
	int i;
	float complex t;

	while (conted == 1) {
		conted = 0;

		switch (ord) {
		case 1:
			for (i = 0; i < N-1; i++) {
				if (GT_(v[i*inc], v[i*inc+inc])) {
					conted = 1;
					SWAP_(v[i*inc], v[i*inc+inc], t);
				}
			}
			break;
		case -1:
			for (i = 0; i < N-1; i++) {
				if (GT_(v[i*inc+inc], v[i*inc])) {
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
void cmerge_(
	int N,
	int N1, float complex* v1, int inc1, int* ptr_id1,
	int N2, float complex* v2, int inc2, int* ptr_id2,
	float complex* v, int inc,
	int ord
) {
	int id1, id2, id;

	id1 = ptr_id1[0];
	id2 = ptr_id2[0];

	switch (ord) {
		case 1:
			for (id = 0; id < N && id1 < N1 && id2 < N2; id++) {
				if (GT_(v1[id1*inc1], v2[id2*inc2])) {
					v[id*inc] = v1[inc1*(id1++)];
				}
				else
					v[id*inc] = v2[inc2*(id2++)];
			}
			break;
		case -1:
			for (id = 0; id < N && id1 < N1 && id2 < N2; id++) {
				if (GT_(v2[id1*inc], v1[id2*inc])) {
					v[id*inc] = v1[inc1*(id1++)];
				}
				else
					v1[id*inc] = v2[inc2*(id2++)];
			}
			break;
		case 2:
			for (id = 0; id < N && id1 < N1 && id2 < N2; id++) {
				if (fabsz(v1[id1*inc1]) < fabsz(v2[id2*inc2])) {
					v[id*inc] = v1[inc1*(id1++)];
				}
				else
					v[id*inc] = v2[inc2*(id2++)];
			}
			break;
		case -2:
			for (id = 0; id < N && id1 < N1 && id2 < N2; id++) {
				if (fabsz(v1[id1*inc1]) > fabsz(v2[id2*inc2])) {
					v[id*inc] = v1[inc1*(id1++)];
				}
				else
					v[id*inc] = v2[inc2*(id2++)];
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

void cmerge_sort_(
	int N1, float complex* v1, int inc1,
	int N2, float complex* v2, int inc2,
	int ord,
	float complex* buffer
) {
	int id1, id2;
	int i;

	id1 = id2 = 0;
	//if (N1 <= N2) {
		// COPY v1 to buffer
		for (i = 0; i < N1; i++)
			buffer[i] = v1[i * inc1];

		// FIRST PHASE: MERGE BUFFER,V2 TO V1
		cmerge_(N1, N1, buffer, 1, &id1, N2, v2, inc2, &id2, v1, inc1, ord);

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
		cmerge_(N2, N1, buffer, 1, &id1, N2, v2, inc2, &id2, v2, inc2, ord);

		return;
	//}
}

void csort_merge_(
	int N,
	float complex* v, int inc,
	int ord,
	float complex* buffer
) {
	if (N <= 32) {
		csort_bubble(N, v, inc, ord);
		return;
	}

	int mid = N /2;

	// left recursive call
	csort_merge_(mid, v, inc, ord, buffer);
	// right recursive call
	csort_merge_(N-mid, v+mid*inc, inc, ord, buffer);
	// merge
	cmerge_sort_(mid, v, inc, N-mid, v+mid*inc, inc, ord, buffer);
}

void csort_merge(
	int N,
	float complex* v, int inc,
	int ord
) {
	float complex* buffer = (float complex*) malloc((N/2) * sizeof(float complex));

	csort_merge_(N, v, inc, ord, buffer);

	free(buffer);
}

