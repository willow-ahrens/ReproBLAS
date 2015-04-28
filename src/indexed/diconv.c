/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

void dmdconv(const int fold, double x, double* repy, int increpy, double* cary, int inccary) {
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
	dmbound(fold, dindex(fabs(x)), repy, increpy);
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

void didconv(const int fold, double x, double_indexed *y) {
	dmdconv(fold, x, y, 1, y + fold, 1);
}

void zmzconv(const int fold, void *x, double *repy, int increpy, double *cary, int inccary) {
  dmdconv(fold, ((double*)x)[0], repy, increpy * 2, cary, inccary * 2);
  dmdconv(fold, ((double*)x)[1], repy + 1, increpy * 2, cary + 1, inccary * 2);
}

void zizconv(const int fold, void *x, double_complex_indexed *y) {
	zmzconv(fold, x, y, 1, y + 2 * fold, 1);
}

double ddmconv(const int fold, double* repx, int increpx, double* carx, int inccarx) {
	int i;
	double y = 0.0;

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

double ddiconv(const int fold, double_indexed *x) {
  return ddmconv(fold, x, 1, x + fold, 1);
}

void zzmconv_sub(const int fold, double *repx, int increpx, double *carx, int inccarx, void *y) {
	((double*)y)[0] = ddmconv(fold, repx, increpx * 2, carx, inccarx + 1);
	((double*)y)[1] = ddmconv(fold, repx + 1, increpx * 2, carx + 1, inccarx + 1);
}

void zziconv_sub(const int fold, double_complex_indexed *x, void *y) {
  zzmconv_sub(fold, x, 1, x + 2 * fold, 1, y);
}
