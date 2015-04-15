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

#define CHECK_NAN_INF
#define POW2(x) powf(2.0,x)

// PRE-FIXED BOUNDARIES
#define SDEFAULT_W 15
#define PREC 24 
#define F_MMASK ~((1 << PREC) - 1)
#define F_BOUDNARY_ZERO_IND 10
static float F_BOUNDARIES[20];
static int   F_BOUNDARIES_initialized = 0;
static int   F_BOUNDARY_NB    = 1 << (PREC - 2 - SDEFAULT_W);
static int   F_BIN_WIDTH      = SDEFAULT_W;
static float F_BOUNDARY_STEP  = (float)(1 << SDEFAULT_W);
static float F_BOUNDARY_STEP1 = 1.0/(1 << SDEFAULT_W);
static int   F_BOUNDARY_LEFT  = 10;
static int   F_BOUNDARY_RIGHT = 10;

#define FF_ 1.5

int sIWidth() { return F_BIN_WIDTH; }
int sICapacity() { return F_BOUNDARY_NB; }

void sI_Initialize_() {
	if (F_BOUNDARIES_initialized) return;
	F_BOUNDARIES[F_BOUDNARY_ZERO_IND] = FF_;
	int exp = -1;
	int ind = F_BOUDNARY_ZERO_IND + 1;
	float bound = F_BOUNDARY_STEP1 * FF_;
	while (exp * F_BIN_WIDTH  >= FLT_MIN_EXP) {
		F_BOUNDARIES[ind] = bound;
		ind++;
		exp--;
		bound *= F_BOUNDARY_STEP1;
	}
	F_BOUNDARY_RIGHT = ind;
	while (ind < 20) {
		F_BOUNDARIES[ind] = 0.0;
		ind++;
	}

	exp = 1;
	bound = F_BOUNDARY_STEP * FF_;
	ind = F_BOUDNARY_ZERO_IND - 1;
	while (exp * F_BIN_WIDTH <= FLT_MAX_EXP) {
		F_BOUNDARIES[ind] = bound;
		ind--;
		exp++;
		bound *= F_BOUNDARY_STEP;
	}
	F_BOUNDARY_LEFT = ind;
	F_BOUNDARIES[ind--] = bound;

	while (ind >= 0) {
		F_BOUNDARIES[ind] = INFINITY;
		ind--;
	}
	F_BOUNDARIES_initialized = 1;
}

float F_Ind2Boundary(int index) {
	sI_Initialize_();
	index = F_BOUDNARY_ZERO_IND - index;
	return F_BOUNDARIES[index];
}

float* F_Ind2Boundares(int index) {
	sI_Initialize_();
	index = F_BOUDNARY_ZERO_IND - index;
	return F_BOUNDARIES+index;
}

int F_Max2Ind(float amax) {
	sI_Initialize_();
	if (amax == 0)
		return F_BOUNDARY_RIGHT;
	if (isinf(amax))
		return F_BOUNDARY_LEFT;
	
	int mid = F_BOUDNARY_ZERO_IND;
	
	while (F_BOUNDARIES[mid] > amax) mid++;
	mid--;
	while (F_BOUNDARIES[mid] < amax) mid--;

	// TODO: CHECK ...
	if (amax * (4 * FF_) * F_BOUNDARY_NB > F_BOUNDARIES[mid])
		mid--;
	return mid;
}

float* F_Max2Boundaries(float amax) {
	return F_BOUNDARIES + F_Max2Ind(amax);
}

// COMPUTE THE BOUNDARIES BASED ON MAXIMUM ABSOLUTE VALUE
int sIBoundary_(int fold, float max, float* M, int inc) {
	float delta;
	int i;
	float M0;

	int index;
    index = F_Max2Ind(max);
    for (i = 0; i < fold; i++, M += inc) {
        M[0] = F_BOUNDARIES[index + i];
    }
    return index;

    /*
	if (step == F_BIN_WIDTH || step == 0) {
		index = F_Max2Ind(max);
		for (i = 0; i < fold; i++, M += inc) {
			M[0] = F_BOUNDARIES[index + i];
		}
		return index;
	}

	// TODO: CHECK FOR INFINITY & NAN
	int log2N = PREC - 2 - step;
	delta = ceil(log2(max));
	index = (int) (delta) + log2N + 2;
	index = ((step * 1024 + index + step - 1) / step - 1024);

	index *= step;
	float dstep = powf(2.0, -step);

	if (index >= FLT_MAX_EXP) {
		// TO AVOID FALSE OVERFLOW
		M0 = POW2( index - step);
		M0 *= FF_;
		M[0] = M0;
		if (fold > 1)
			M[1] = M0;
	}
	else {
		M0 = POW2( index);
		M0 *= FF_;
		M[0] = M0;
		if (fold > 1) {
			M0 *= dstep;
			M[inc] = M0;
		}
	}

	for (i = 2; i < fold; i++) {
		M0   *=  dstep;
		M[i*inc] = M0;
	}

#ifdef DEBUG
	fprintf(stdout, "EXTRACTION FACTORS: [%g", M[0]);
	for (i = 1; i < fold; i++)
		fprintf(stdout, ",%g", M[i*inc]);
	fprintf(stdout, "]\n");
#endif

	return index / step;
    */
}

#define USE_FREXP
float ufpf(float x) {
	if (x == 0.0)
		return 0.0;
#ifdef USE_FREXP
	int exp;
	frexpf(x, &exp);
	return ldexpf(0.5, exp);
#elif defined ( __SSE__ )
	__m128 mX, mMask;
	mX = _mm_load_ss(&x);
	mMask = _mm_set_ss(INFINITY);
	mX = _mm_and_ps(mX, mMask);
	_mm_store_ss(&x, mX);
	return x;
#else
	int_float lM;
	lM.f = x;
	lM.i &= F_MMASK;
	return lM.f;
#endif
}

void sIprint1(int n, float* x, float* carry, int inc) {
	int i;
	float M;
	for (i = 0; i < n; i++) {
		M = ufpf(x[i*inc]);
		printf("{M:2^%g # %g :: %g (%g)}", log2f(M), carry[i*inc], x[i*inc] - 1.5*M, x[i*inc] + (carry[i*inc]-6)*0.25*M);
	}
}

void cIprint1(int n, float complex* x, float* carry, int inc) {
	/*
	printf("\nReal: ");
	sIprint_(n, (float*) x, carry, inc*2);
	printf("\nImaginary: ");
	sIprint_(n, ((float*) x) + 1, carry+1, inc*2);
	*/
	int i;
	float M;
	float* ptr = (float*) x;
	for (i = 0; i < n; i++, ptr += 2 * inc, carry += 2 * inc) {
		M = ufpf(ptr[0]);
		printf("M:2^%2g", log2f(M));
		printf("# %4g", carry[0]);
		printf(" %8.3g (%8.3g) ", ptr[0], ptr[0] - 1.5 * M);

		M = ufpf(ptr[1]);
		printf(" || M:2^%2g", log2f(M));
		printf("# %6g", carry[1]);
		printf(" %8.3g (%8.3g) ", ptr[1], ptr[1] - 1.5 * M);
		printf("\n");
	}
}

