#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "debug.h"

#define SWAP_(A,B,T) T = A; A = B; B = T;

void drandomseed() {
	struct timeval st;
      	gettimeofday( &st, NULL );
	srand48((long)(st.tv_usec + 1e6*st.tv_sec));
}

double dsum(int N, double* v) {
	double ret = 0.0;
	int i;
	for (i = 0; i < N; i++)
		ret += v[i];
	return ret;
}

void dreverse_(int N, double* v, int inc) {
	int i;
	double t;
	for (i = 0 ; i < N / 2; i++)
	{
		SWAP_(v[i*inc], v[(N -i-1)*inc], t);	
	}
}
void dreverse(int N, double* v) {
	int i;
	double t;
	for (i = 0 ; i < N / 2; i++)
	{
		SWAP_(v[i], v[N -i-1], t);	
	}
}

void dshuffle_(int N, double* v, int inc) {
	int i;
	double t;
	int gap;
	for (i = 0 ; i < N - 2; i++)
	{
		gap = rand() % (N - i);
		SWAP_(v[i*inc], v[(i + gap)*inc], t);	
	}
}
void dshuffle(int N, double* v) {
	int i;
	double t;
	int gap;
	for (i = 0 ; i < N - 2; i++)
	{
		gap = rand() % (N - i);
		SWAP_(v[i], v[i + gap], t);	
	}
}

double pdsum(int N, double* v, int procs) {
	int lN = (N + procs - 1) / procs;
	int i;
	double ret = 0.0;
	int n;

	for (i = 0; i < N; i += lN, v+=lN) {
		n = lN > (N - i) ? (N-i):lN;
		ret += dsum(n, v);
	}
	return ret;
}

void drandcond_(int N, double* v, int inc, double cond) {
	int quart = (N / 8) & ~1;
	int mid = N - quart * 2;
	double c1, c2, c, f;
	int i;

	dgenvec_(quart / 2, v, 0, inc, 1.0e-10);
	dgenvec_(quart / 2, v + inc*(quart / 2), inc, 0, 1);
	dgenvec_(quart / 8, v, inc, 0, 1e-20);
	c1 = 0.0;
	for (i = 0; i < quart; i++) {
		c1 += v[i*inc];
		v[(i + quart)*inc] = -v[i*inc];
	}
	dgenvec_(quart, v + 2 * quart*inc, inc, 0, 1.0);
	dgenvec_(mid - quart, v + 3 * quart *inc, inc, 0, 1e-8);
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

	dshuffle_(N, v, inc);
}
void drandcond(int N, double* v, double cond) {
	drandcond_(N, v, 1, cond);
}

#define RAND (drand48() * (1+1e-9))
void dgenvec_(int N, double* v, int inc, int type, double factor) {
	int i;
	double v1;
	if (type == 0) {
		for (i = 0; i < N; i++) v[i*inc] = factor * RAND;
	}
		
	if (type == 1) {
		for (i = 0; i < N; i++) v[i*inc] = factor * (2 * RAND - 1);
	}
	if (type == 2) {
		for (i = 0; i < N; i++) v[i*inc] = RAND ;
		for (i = 0; i < N; i++) {
			v1 = RAND - 1;
			v[i*inc] += v1;
		}
	}
	if (type == 3) {
		for (i = 0; i < N; i++) {
			double t1 = drand48();
			double t2 = drand48();
			v1 = sqrt(-2.0 * log(t1)) * cos(2.0 * M_PI * t2);
			v[i * inc] = v1;
		}
	}
	if (type == 4) {
		for (i = 0; i < N; i++) {
			v1 = factor * sin(2.0 * M_PI * ((double)i / (double)N));
			v[i*inc] = v1;
		}
	}
	if (type == 5) {
		drandcond(N, v, v[0]);
	}

	if (type == 10 || type == 11 || type == 12 || type == 13){
		dgenvec_(N/2, v, inc, type - 10, factor);
		dgenvec_(N/2, v+(N/2) * inc, inc, type - 10, factor * 1e-12);
	}

	if (factor != 1.0) {
		for (i = 0; i < N; i++) v[i*inc] *= factor;
	}
}

void dgenvec(int N, double* v, int type, double factor) {
	dgenvec_(N, v, 1, type, factor);
}
// SLOW SORTING
//    ord: 1  INCREASING
//         -1 DECREASING
//         2  INCREASING MAGNITUDE
//         -2 DECREASING MAGNITUDE
void dsort_bubble(int N, double* v, int ord) {
	int conted = 1;
	int i;
	double t;

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
void dmerge_(
	int N,
	int N1, double* v1, int* ptr_id1,
	int N2, double* v2, int* ptr_id2,
	double* v,
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
			memcpy(v + id, v2 + id2, (N - id) * sizeof(double));
		id2 += N - id;
	}
	else if (id2 == N2 && id < N) {
		if (v != v1)
			memcpy(v + id, v1 + id1, (N - id) * sizeof(double));
		id1 += N - id;
	}

	ptr_id1[0] = id1;
	ptr_id2[0] = id2;
}

void dmerge_sort_(
	int N1, double* v1,
	int N2, double* v2,
	int ord,
	double* buffer
) {
	int id1, id2;

	id1 = id2 = 0;
	//if (N1 <= N2) {
		// COPY v1 to buffer
		memcpy(buffer, v1, N1 * sizeof(double));

		// FIRST PHASE: MERGE BUFFER,V2 TO V1
		dmerge_(N1, N1, buffer, &id1, N2, v2, &id2, v1, ord);

		// buffer empty: done
		if (id1 == N1)
			return;

		// V2 empty: copy buffer to v2
		if (id2 == N2) {
			memcpy(v2, buffer, N2 * sizeof(double));
			return;
		}

		//SECOND PHASE: MERGE BUFFER+ID1, V2+ID2 TO V2
		dmerge_(N2, N1, buffer, &id1, N2, v2, &id2, v2, ord);

		return;
	//}
}

void dsort_merge_(
	int N,
	double* v,
	int ord,
	double* buffer
) {
	if (N <= 32) {
		dsort_bubble(N, v, ord);
		return;
	}

	int mid = N /2;

	// left recursive call
	dsort_merge_(mid, v, ord, buffer);
	// right recursive call
	dsort_merge_(N-mid, v+mid, ord, buffer);
	// merge
	dmerge_sort_(mid, v, N-mid, v+mid, ord, buffer);
}

void dsort_merge(
	int N,
	double* v,
	int ord
) {
	double* buffer = (double*) malloc((N/2) * sizeof(double));

	dsort_merge_(N, v, ord, buffer);

	free(buffer);
}

