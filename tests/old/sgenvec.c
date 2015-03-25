#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "debug.h"

#define SWAP_(A,B,T) T = A; A = B; B = T;

float ssum(int N, float* v) {
	float ret = 0.0;
	int i;
	for (i = 0; i < N; i++)
		ret += v[i];
	return ret;
}

void sreverse_(int N, float* v, int inc) {
	int i;
	float t;
	for (i = 0 ; i < N / 2; i++)
	{
		SWAP_(v[i*inc], v[(N -i-1)*inc], t);	
	}
}
void sreverse(int N, float* v) {
	sreverse_(N,v,1);
}

void sshuffle_(int N, float* v, int inc) {
	int i;
	float t;
	int gap;
	for (i = 0 ; i < N - 2; i++)
	{
		gap = rand() % (N - i);
		SWAP_(v[i*inc], v[(i + gap)*inc], t);	
	}
}
void sshuffle(int N, float* v) {
	sshuffle_(N,v,1);
}

float pssum(int N, float* v, int procs) {
	int lN = (N + procs - 1) / procs;
	int i;
	float ret = 0.0;
	int n;

	for (i = 0; i < N; i += lN, v+=lN) {
		n = lN > (N - i) ? (N-i):lN;
		ret += ssum(n, v);
	}
	return ret;
}

void srandcond(int N, float* v, int inc, float cond) {
	int quart = (N / 8) & ~1;
	int mid = N - quart * 2;
	double c1, c2, c, f;
	int i;

	sgenvec_(quart / 2, v, inc, 0, 1.0e-10);
	sgenvec_(quart / 2, v + (quart / 2) * inc, inc, 0, 1);
	sgenvec_(quart / 8, v, inc, 0, 1e-20);
	c1 = 0.0;
	for (i = 0; i < quart; i++) {
		c1 += v[i*inc];
		v[i + quart] = -v[i*inc];
	}
	sgenvec_(quart, v + 2 * quart * inc, inc, 0, 1.0);
	sgenvec_(mid - quart, v + 3 * quart * inc, inc, 0, 1e-8);
	c2 = 0.0; c = 0.0;
	for (i = 2 * quart; i < N; i++) {
		c2 += fabs(v[i*inc]);
		c  += v[i*inc];
	}

	f = 2 * c1 / (cond * fabs(c) - c2);
	if (f < 1e-16) f = 1e-16;
	if (f > 1e16) f = 1e16;
	for (i = 2 * quart; i < N; i++) {
		v[i*inc] *= f;
	}

	sshuffle_(N, v, 1);
}

#define RAND ((float)drand48() * (1+1e-4))
void sgenvec_(int N, float* v, int inc, int type, float factor) {
	int i;
	double t1, t2, v1;
	if (type == 0) {
		for (i = 0; i < N; i++) v[i*inc] = factor * RAND;
	}

	if (type == 1) {
		for (i = 0; i < N; i++) v[i*inc] = factor * (2 * RAND - 1);
	}
	if (type == 2) {
		for (i = 0; i < N; i++) v[i*inc] = RAND ;
		for (i = 0; i < N; i++) {
			v[i*inc] += RAND - 1;
		}
	}
	if (type == 3) {
		for (i = 0; i < N; i++) {
			t1 = drand48();
			t2 = drand48();
			v1 = sqrt(-2.0 * log(t1)) * cos(2.0 * M_PI * t2);
			v[i*inc] = (float)v1;
		}
	}
	if (type == 4) {
		for (i = 0; i < N; i++) {
			v1 = factor * sin(2.0 * M_PI * ((float)i / (float)N));
			v[i * inc] = (float) v1;
		}
	}
	if (type == 5) {
		srandcond(N, v, inc, v[0]);
	}

	if (type == 10 || type == 11 || type == 12 || type == 13){
		sgenvec_(N/2, v, inc, type - 10, factor);
		sgenvec_(N/2, v+(N/2) * inc, inc, type - 10, factor * 1e-12);
	}

	if (factor != 1.0) {
		for (i = 0; i < N; i++) v[i*inc] *= factor;
	}
}
void sgenvec(int N, float* v, int type, float factor) {
	sgenvec_(N, v, 1, type, factor);
}

// SLOW SORTING
//    ord: 1  INCREASING
//         -1 DECREASING
//         2  INCREASING MAGNITUDE
//         -2 DECREASING MAGNITUDE
void ssort_bubble(int N, float* v, int ord) {
	int conted = 1;
	int i;
	float t;

	while (conted == 1) {
		conted = 0;

		switch (ord) {
		case 1:
			for (i = 0; i < N-1; i++) {
				if (v[i] > v[i+1]) {
					conted = 1;
					SWAP_(v[i], v[i+1], t);
				}
			}
			break;
		case -1:
			for (i = 0; i < N-1; i++) {
				if (v[i] < v[i+1]) {
					conted = 1;
					SWAP_(v[i], v[i+1], t);
				}
			}
			break;
		case 2:
			for (i = 0; i < N-1; i++) {
				if (fabs(v[i]) > fabs(v[i+1])) {
					conted = 1;
					SWAP_(v[i], v[i+1], t);
				}
			}
			break;
		case -2:
			for (i = 0; i < N-1; i++) {
				if (fabs(v[i]) < fabs(v[i+1])) {
					conted = 1;
					SWAP_(v[i], v[i+1], t);
				}
			}
			break;
		}
	}
}

// merge v1, v2 into v
void smerge_(
	int N,
	int N1, float* v1, int* ptr_id1,
	int N2, float* v2, int* ptr_id2,
	float* v,
	int ord
) {
	int id1, id2, id;

	id1 = ptr_id1[0];
	id2 = ptr_id2[0];

	switch (ord) {
		case 1:
			for (id = 0; id < N && id1 < N1 && id2 < N2; id++) {
				if (v1[id1] < v2[id2]) {
					v[id] = v1[id1++];
				}
				else
					v[id] = v2[id2++];
			}
			break;
		case -1:
			for (id = 0; id < N && id1 < N1 && id2 < N2; id++) {
				if (v1[id1] > v2[id2]) {
					v[id] = v1[id1++];
				}
				else
					v1[id] = v2[id2++];
			}
			break;
		case 2:
			for (id = 0; id < N && id1 < N1 && id2 < N2; id++) {
				if (fabs(v1[id1]) < fabs(v2[id2])) {
					v[id] = v1[id1++];
				}
				else
					v[id] = v2[id2++];
			}
			break;
		case -2:
			for (id = 0; id < N && id1 < N1 && id2 < N2; id++) {
				if (fabs(v1[id1]) > fabs(v2[id2])) {
					v[id] = v1[id1++];
				}
				else
					v[id] = v2[id2++];
			}
			break;
	}

	if (id1 == N1 && id < N) {
		if (v != v2)
			memcpy(v + id, v2 + id2, (N - id) * sizeof(float));
		id2 += N - id;
	}
	else if (id2 == N2 && id < N) {
		if (v != v1)
			memcpy(v + id, v1 + id1, (N - id) * sizeof(float));
		id1 += N - id;
	}

	ptr_id1[0] = id1;
	ptr_id2[0] = id2;
}

void smerge_sort_(
	int N1, float* v1,
	int N2, float* v2,
	int ord,
	float* buffer
) {
	int id1, id2;

	id1 = id2 = 0;
	//if (N1 <= N2) {
		// COPY v1 to buffer
		memcpy(buffer, v1, N1 * sizeof(float));

		// FIRST PHASE: MERGE BUFFER,V2 TO V1
		smerge_(N1, N1, buffer, &id1, N2, v2, &id2, v1, ord);

		// buffer empty: done
		if (id1 == N1)
			return;

		// V2 empty: copy buffer to v2
		if (id2 == N2) {
			memcpy(v2, buffer, N2 * sizeof(float));
			return;
		}

		//SECOND PHASE: MERGE BUFFER+ID1, V2+ID2 TO V2
		smerge_(N2, N1, buffer, &id1, N2, v2, &id2, v2, ord);

		return;
	//}
}

void ssort_merge_(
	int N,
	float* v,
	int ord,
	float* buffer
) {
	if (N <= 32) {
		ssort_bubble(N, v, ord);
		return;
	}

	int mid = N /2;

	// left recursive call
	ssort_merge_(mid, v, ord, buffer);
	// right recursive call
	ssort_merge_(N-mid, v+mid, ord, buffer);
	// merge
	smerge_sort_(mid, v, N-mid, v+mid, ord, buffer);
}

void ssort_merge(
	int N,
	float* v,
	int ord
) {
	float* buffer = (float*) malloc((N/2) * sizeof(float));

	ssort_merge_(N, v, ord, buffer);

	free(buffer);
}

