/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <fenv.h>

#include "dIndexed.h"

// a += b
void ddpdd(double* a, double* b) {
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = a[0] + b[0];
	  bv = s1 - a[0];
	  s2 = ((b[0] - bv) + (a[0] - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = a[1] + b[1];
	  bv = t1 - a[1];
	  t2 = ((b[1] - bv) + (a[1] - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  a[0] = t1 + t2;
	  a[1] = t2 - (a[0] - t1);
}

void ddpd(double* a, double b) {
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = a[0] + b;
	  bv = s1 - a[0];
	  s2 = ((b - bv) + (a[0] - (s1 - bv)));

	  t1 = a[1] + s2;
	  bv = t1 - a[1];
	  t2 = ((s2 - bv) + (a[1] - (t1 - bv)));

	  s2 = t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  t2 += s2 - (t1 - s1);

	  /* Renormalize (t1, t2)  */
	  a[0] = t1 + t2;
	  a[1] = t2 - (a[0] - t1);
}

void dndpd(int fold, double* a, double b) {
	double bv;
	double s1, s2;
	int i, j;

	// DISTILATION
	for (i = 0; i < fold; i++) {
		/* Add two hi words. */
		s1 = a[i] + b;
		bv = s1 - a[i];
		s2 = ((b - bv) + (a[i] - (s1 - bv)));

		a[i] = s1;
		b = s2;
	}
	for (j = 0; j < fold + 1; j++) {
		for (i = 0; i < fold - 1; i++) {
			/* Add two hi words. */
			s1 = a[i] + a[i+1];
			bv = s1 - a[i];
			s2 = ((a[i+1] - bv) + (a[i] - (s1 - bv)));

			a[i] = s1;
			a[i+1] = s2;
		}
		/* Add two hi words. */
		s1 = a[fold - 1] + b;
		bv = s1 - a[fold - 1];
		s2 = ((b - bv) + (a[fold - 1] - (s1 - bv)));

		a[fold - 1] = s1;
		b = s2;
	}
}

#define SCALE_ 67108864 // 2^26
double TwoProd(double a, double b, double* lo) {
	double r = a * b;
	double tmp, ah, al, bh, bl;

	if (lo == NULL)
		return r;

	tmp = a * SCALE_;
	ah = tmp + a;
	ah = ah - tmp;
	al = a - ah;

	tmp = b * SCALE_;
	bh = tmp + b;
	bh = bh - tmp;
	bl = b - bh;

	*lo = (((ah * bh - r) + ah * bl) + al * bh) + al * bl;

	return r;
}

double TwoSum(double a, double b, double* lo) {
	double r = a + b;
	double s;

	if (lo == NULL)
		return r;

	s = r - a;
	*lo = ((b - s) + (a - (a - s)));

	return r;
}

