/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

void dmdconv(double x, double* repy, int increpy, double* cary, int inccary, int fold) {
	int i;
	double q;
	double s;
	if (x == 0.0) {
		for (i = 0; i < fold; i++) {
			repy[i*increpy] = 0.0;
			cary[i*inccary] = 0.0;
		}
		return;
	}
	dIBoundary(fold, fabs(x), repy, increpy);
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

void didconv(double x, double_indexed *y, int fold) {
	dmdconv(x, y, 1, y + fold, 1, fold);
}

void zmzconv(void *x, double *repy, int increpy, double *cary, int inccary, int fold) {
  dmdconv(((double*)x)[0], repy, increpy * 2, cary, inccary * 2, fold);
  dmdconv(((double*)x)[1], repy + 1, increpy * 2, cary + 1, inccary * 2, fold);
}

void zizconv(void *x, double_complex_indexed *y, int fold) {
	zmzconv(x, y, 1, y + 2 * fold, 1, fold);
}

double ddmconv(double* repx, int increpx, double* carx, int inccarx, int fold) {
	int i;
	double y = 0.0;

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

double ddiconv(double_indexed *x, int fold) {
  return ddmconv(x, 1, x + fold, 1, fold);
}

void zzmconv_sub(double *repx, int increpx, double *carx, int inccarx, void *y, int fold) {
	((double*)y)[0] = ddmconv(repx, increpx * 2, carx, inccarx + 1, fold);
	((double*)y)[1] = ddmconv(repx + 1, increpx * 2, carx + 1, inccarx + 1, fold);
}

void zziconv_sub(double_complex_indexed *x, void *y, int fold) {
  zzmconv_sub(x, 1, x + 2 * fold, 1, y, fold);
}
