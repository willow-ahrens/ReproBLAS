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
#define WIDTH 64
#define BIN_WIDTH 40
#define PREC 53
#define BOUND_ZERO_INDEX 32
static double bounds[64];
static int    bounds_initialized = 0;
static int    capacity    = 1 << (PREC - BIN_WIDTH - 2); // 2^(51-W)
static double bound_step;//initialized in bounds_intialize
static int    bound_min_index  = 32;//initialized in bounds_intialize
static int    bound_max_index = 32;//initialized in bounds_intialize

int dIWidth() { return BIN_WIDTH; }
int dICapacity() { return capacity; }

static void bounds_initialize() {
	if (bounds_initialized) return;
    bound_step = ldexp(1, BIN_WIDTH);
	bounds[BOUND_ZERO_INDEX] = 1.5;
	int exp = -1;
	int ind = BOUND_ZERO_INDEX + 1;
	double bound = (1.0/bound_step) * 1.5;
	while (exp * BIN_WIDTH  >= DBL_MIN_EXP) {
		bounds[ind] = bound;
		ind++;
		exp--;
		bound /= bound_step;
	}
	bound_max_index = ind;
	while (ind < 64) {
		bounds[ind] = 0.0;
		ind++;
	}

	exp = 1;
	bound = bound_step * 1.5;
	ind = BOUND_ZERO_INDEX - 1;
	while (exp * BIN_WIDTH <= DBL_MAX_EXP) {
		bounds[ind] = bound;
		ind--;
		exp++;
		bound *= bound_step;
	}
	bound_min_index = ind;
	while (ind >= 0) {
		bounds[ind] = bound;
		ind--;
	}
	bounds_initialized = 1;
}

/*
double D_Ind2Boundary(int index) {
	bounds_initialize();
	index = BOUND_ZERO_INDEX - ((index & 2047) - 1024);
	return bounds[index];
}

double* D_Ind2Boundaries(int index) {
	bounds_initialize();
	index = BOUND_ZERO_INDEX - ((index & 2047) - 1024);
	return bounds+index;
}
*/

int D_Max2Ind(double amax) {
	bounds_initialize();
	if (amax == 0)
		return bound_max_index;
	// TODO OVERFLOW
	amax *= 1.5 * (1 << (PREC - BIN_WIDTH));
	int left = bound_min_index;
	int right = bound_max_index;
	int mid = BOUND_ZERO_INDEX;

	while (right - left > 0) {
		// FOUND
		if (bounds[mid+1] <= amax && amax < bounds[mid]) {
			return mid;
		}
		
		if (bounds[mid] <= amax) {
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
	bounds_initialize();
	if (amax == 0)
		return bounds + bound_max_index;
	// TODO OVERFLOW
	amax *= 1.5 * capacity;
	int left = bound_min_index;
	int right = bound_max_index;
	int mid = BOUND_ZERO_INDEX;

	while (right - left > 0) {
		// FOUND
		if (bounds[mid+1] < amax && amax <= bounds[mid]) {
			return bounds + mid;
		}
		
		if (bounds[mid] < amax) {
			right = mid;
		}
		else {
			left = mid;
		}
		mid = (left + right) / 2;
	}
	return bounds + left;
}


// COMPUTE THE BOUNDARIES BASED ON MAXIMUM ABSOLUTE VALUE
int dIBoundary(int fold, double max, double* M, int inc) {
	int i;
    int index;
    int other = D_Max2Ind(max);

    max *= 1.5 * (1l << (PREC - BIN_WIDTH));
    if(max >= bounds[bound_min_index]){
      index = bound_min_index;
    } else if(max < bounds[bound_max_index]){
      index = bound_max_index;
    } else {
      if(frexp(max, &index) < 0.75){
        index--;
      }
      index--;
      if(index < 0){
        index -= BIN_WIDTH - 1; //fixing integer division rounding negatives towards 0
      }
      index /= BIN_WIDTH;
      index = BOUND_ZERO_INDEX - 1 - index;
    }
    if(index != other){
      printf("max %g, index %d, peter index %d %g\n", max, other, index, bounds[index]);
    }
    for (i = 0; i < fold; i++) {
        M[i * inc] = bounds[index + i];
    }
    return index;


    /*
	if (W == BIN_WIDTH || W == 0) {
		index = D_Max2Ind(max);
		for (i = 0; i < fold; i++) {
			M[i * inc] = bounds[index + i];
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
		M0 *= 1.5;
		M[0] = M0;
		if (fold > 1)
			M[inc] = M0;
	}
	else {
		M0 = ldexp(0.5, 1 + index);
		M0 *= 1.5;
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
  */
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

