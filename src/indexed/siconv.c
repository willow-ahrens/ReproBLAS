#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

void smsconv(const int fold, float x, float* repy, int increpy, float* cary, int inccary) {
	int i;
	float q;
	float s;
	if (x == 0.0) {
		for (i = 0; i < fold; i++) {
			repy[i*increpy] = 0.0;
			cary[i*inccary] = 0.0;
		}
		return;
	}
	smbound(fold, sindex(fabs(x)), repy, increpy);
	for (i = 0; i < fold; i++, repy += increpy, cary += inccary) {
		// high-order part
		cary[0] = 0.0;

		s = repy[0];
		q = s + x;
		repy[0] = s;
		q -= s;
		x -= q;
	}
}

void sisconv(const int fold, float x, float_indexed *y) {
	smsconv(fold, x, y, 1, y + fold, 1);
}

void cmcconv(const int fold, void *x, float *repy, int increpy, float *cary, int inccary) {
  smsconv(fold, ((float*)x)[0], repy, increpy * 2, cary, inccary * 2);
  smsconv(fold, ((float*)x)[1], repy + 1, increpy * 2, cary + 1, inccary * 2);
}

void cicconv(const int fold, void *x, float_complex_indexed *y) {
	cmcconv(fold, x, y, 1, y + 2 * fold, 1);
}

float ssmconv(const int fold, float* repx, int increpx, float* carx, int inccarx) {
	int i;
	float y = 0.0;

	// CHECK FOR NAN OR INFINITY
	if (isinf(repx[0]) || isnan(repx[0]))
		return repx[0];

	if (repx[0] == 0.0) {
		return 0.0;
	}

	// TODO: SCALING TO AVOID OVERFLOW

	for (i = 0; i < fold; i++, repx += increpx, carx += inccarx) {
		y += (repx[0] + (carx[0] - 6) * ufp(repx[0]) * 0.25);
	}

	return y;
}

float ssiconv(const int fold, float_indexed *x) {
  return ssmconv(fold, x, 1, x + fold, 1);
}

void ccmconv_sub(const int fold, float *repx, int increpx, float *carx, int inccarx, void *y) {
	((float*)y)[0] = ssmconv(fold, repx, increpx * 2, carx, inccarx + 1);
	((float*)y)[1] = ssmconv(fold, repx + 1, increpx * 2, carx + 1, inccarx + 1);
}

void cciconv_sub(const int fold, float_complex_indexed *x, void *y) {
  ccmconv_sub(fold, x, 1, x + 2 * fold, 1, y);
}
