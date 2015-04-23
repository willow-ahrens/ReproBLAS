/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

#ifdef __SSE2__
#	include <emmintrin.h>
#endif

#define CHECK_NAN_INF

// PRE-FIXED BOUNDARIES
#define DEFAULT_W 40
#define PREC 53
#define D_BOUNDARY_ZERO_IND 32
static double D_BOUNDARIES[64];
static int    D_BOUNDARIES_initialized = 0;
static int    D_BOUNDARY_NB    = 1<<(PREC - DEFAULT_W - 2); // 2^(51-W)
//static int D_BOUNDARY_LOG2NB = 10;
static int    D_BIN_WIDTH      = DEFAULT_W;
static double D_BOUNDARY_STEP  = (double)(1l << DEFAULT_W);
static double D_BOUNDARY_STEP1 = 1.0/(1l << DEFAULT_W);
static int    D_BOUNDARY_LEFT  = 32;
static int    D_BOUNDARY_RIGHT = 32;

#define FF_ 1.5

int dIWidth() { return D_BIN_WIDTH; }
int dICapacity() { return D_BOUNDARY_NB; }

void dI_initialize_() {
	if (D_BOUNDARIES_initialized) return;
	D_BOUNDARIES[D_BOUNDARY_ZERO_IND] = FF_;
	int exp = -1;
	int ind = D_BOUNDARY_ZERO_IND + 1;
	double bound = D_BOUNDARY_STEP1 * FF_;
	while (exp * D_BIN_WIDTH  >= DBL_MIN_EXP) {
		D_BOUNDARIES[ind] = bound;
		ind++;
		exp--;
		bound *= D_BOUNDARY_STEP1;
	}
	D_BOUNDARY_RIGHT = ind;
	while (ind < 64) {
		D_BOUNDARIES[ind] = 0.0;
		ind++;
	}

	exp = 1;
	bound = D_BOUNDARY_STEP * FF_;
	ind = D_BOUNDARY_ZERO_IND - 1;
	while (exp * D_BIN_WIDTH <= DBL_MAX_EXP) {
		D_BOUNDARIES[ind] = bound;
		ind--;
		exp++;
		bound *= D_BOUNDARY_STEP;
	}
	D_BOUNDARY_LEFT = ind;
	while (ind >= 0) {
		D_BOUNDARIES[ind] = bound;
		ind--;
	}
	D_BOUNDARIES_initialized = 1;
}

double D_Ind2Boundary(int index) {
	dI_initialize_();
	index = D_BOUNDARY_ZERO_IND - ((index & 2047) - 1024);
	return D_BOUNDARIES[index];
}

double* D_Ind2Boundaries(int index) {
	dI_initialize_();
	index = D_BOUNDARY_ZERO_IND - ((index & 2047) - 1024);
	return D_BOUNDARIES+index;
}

int D_Max2Ind(double amax) {
	dI_initialize_();
	if (amax == 0)
		return D_BOUNDARY_RIGHT;
	// TODO OVERFLOW
	amax *= FF_ * (1 << (PREC - D_BIN_WIDTH));
	int left = D_BOUNDARY_LEFT;
	int right = D_BOUNDARY_RIGHT;
	int mid = D_BOUNDARY_ZERO_IND;

	while (right - left > 0) {
		// FOUND
		if (D_BOUNDARIES[mid+1] <= amax && amax < D_BOUNDARIES[mid]) {
			return mid;
		}
		
		if (D_BOUNDARIES[mid] <= amax) {
			right = mid;
		}
		else {
			left = mid;
		}
		mid = (left + right) / 2;
	}
	return left;
}
double* D_Max2Boundaries(double amax) {
	dI_initialize_();
	if (amax == 0)
		return D_BOUNDARIES + D_BOUNDARY_RIGHT;
	// TODO OVERFLOW
	amax *= FF_ * D_BOUNDARY_NB;
	int left = D_BOUNDARY_LEFT;
	int right = D_BOUNDARY_RIGHT;
	int mid = D_BOUNDARY_ZERO_IND;

	while (right - left > 0) {
		// FOUND
		if (D_BOUNDARIES[mid+1] < amax && amax <= D_BOUNDARIES[mid]) {
			return D_BOUNDARIES + mid;
		}
		
		if (D_BOUNDARIES[mid] < amax) {
			right = mid;
		}
		else {
			left = mid;
		}
		mid = (left + right) / 2;
	}
	return D_BOUNDARIES + left;
}


// COMPUTE THE BOUNDARIES BASED ON MAXIMUM ABSOLUTE VALUE
int dIBoundary(int fold, int W, double max, double* M, int inc) {
	double delta;
	int i;
	double M0;

	int index;

	if (W == D_BIN_WIDTH || W == 0) {
		index = D_Max2Ind(max);
		for (i = 0; i < fold; i++) {
			M[i * inc] = D_BOUNDARIES[index + i];
		}
		return index;
	}

	// TODO: CHECK FOR INFINITY & NAN
	int log2N = PREC - 2 - W;
	delta = ceil(log2(max));
	index = (int) (delta) + log2N + 2;
	index = ((W * 1024 + index + W - 1) / W - 1024);

	index *= W;
	double dW = ldexp(0.5, 1-W);

	if (index >= DBL_MAX_EXP) {
		// TO AVOID FALSE OVERFLOW
		M0 = ldexp(0.5, 1 + index - W);
		M0 *= FF_;
		M[0] = M0;
		if (fold > 1)
			M[inc] = M0;
	}
	else {
		M0 = ldexp(0.5, 1 + index);
		M0 *= FF_;
		M[0] = M0;
		if (fold > 1) {
			M0 *= dW;
			M[inc] = M0;
		}
	}

	for (i = 2; i < fold; i++) {
		M0   *=  dW;
		M[i*inc] = M0;
	}

#ifdef DEBUG
	fprintf(stdout, "EXTRACTION FACTORS: [%g", M[0]);
	for (i = 1; i < fold; i++)
		fprintf(stdout, ",%g", M[i]);
	fprintf(stdout, "]\n");
#endif

	return index / W;
}

#define USE_FREXP
double ufp(double x) {
	if (x == 0.0)
		return 0.0;
#ifdef USE_FREXP
	int exp;
	frexp(x, &exp);
	return ldexp(0.5, exp);
#else
	long Mask = ~((1l << 52) - 1);
	long_double lM;
	lM.d = x;
	lM.l &= Mask;
	return lM.d;
#endif
}

void dIprint1(int n, double *x, double *c, int inc) {
	int i;
	double M;
	for (i = 0; i < n; i++, x += inc, c += inc) {
		M = ufp(x[0]);
		printf("{M:2^%g # %g :: %g (%.16g)}", log2(M), c[0], x[0] - 1.5*M,
			(c[0] - 6) * 0.25 * M + x[0]);
	}
}

void zIprint1(int n, double complex* x, double complex* carry, int inc) {
	int i;
	double M;
	double* ptr = (double*) x;
	double* cptr = (double*) carry;
	inc *= 2;
	for (i = 0; i < n; i++, ptr += inc, cptr += inc) {
		M = ufp(ptr[0]);
		printf("M:2^%2g", log2(M));
		printf("# %4g", cptr[0]);
		printf(" %.8g (%8.3g) ", ptr[0], ptr[0] - 1.5 * M);

		M = ufp(ptr[1]);
		printf(" || M:2^%2g", log2(M));
		printf("# %6g", cptr[1]);
		printf(" %.8g (%8.3g) ", ptr[1], ptr[1] - 1.5 * M);
		printf("\n");
	}
}

