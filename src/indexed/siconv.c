#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

void smsconv(float x, float* repy, int increpy, float* cary, int inccary, int fold) {
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
	smbound(sindex(fabs(x)), repy, increpy, fold);
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

void sisconv(float x, float_indexed *y, int fold) {
	smsconv(x, y, 1, y + fold, 1, fold);
}

void cmcconv(void *x, float *repy, int increpy, float *cary, int inccary, int fold) {
  smsconv(((float*)x)[0], repy, increpy * 2, cary, inccary * 2, fold);
  smsconv(((float*)x)[1], repy + 1, increpy * 2, cary + 1, inccary * 2, fold);
}

void cicconv(void *x, float_complex_indexed *y, int fold) {
	cmcconv(x, y, 1, y + 2 * fold, 1, fold);
}

float ssmconv(float* repx, int increpx, float* carx, int inccarx, int fold) {
	int i;
	float y = 0.0;

	// CHECK FOR NAN OR INFINITY
	if (isinf(repx[0]) || isnan(repx[0]))
		return repx[0];

	if (repx[0] == 0.0) {
		return 0.0;
	}

	// TODO: SCALING TO AVOID OVERFLOW
	fold = (fold < 1) ? 1 : fold; //TODO wtf is this line for?

	for (i = 0; i < fold; i++, repx += increpx, carx += inccarx) {
		y += (repx[0] + (carx[0] - 6) * ufp(repx[0]) * 0.25);
	}

	return y;
}

float ssiconv(float_indexed *x, int fold) {
  return ssmconv(x, 1, x + fold, 1, fold);
}

void ccmconv_sub(float *repx, int increpx, float *carx, int inccarx, void *y, int fold) {
	((float*)y)[0] = ssmconv(repx, increpx * 2, carx, inccarx + 1, fold);
	((float*)y)[1] = ssmconv(repx + 1, increpx * 2, carx + 1, inccarx + 1, fold);
}

void cciconv_sub(float_complex_indexed *x, void *y, int fold) {
  ccmconv_sub(x, 1, x + 2 * fold, 1, y, fold);
}
