#include <rblas.h>
#include <IndexedFP.h>
#include <stdio.h>
#include <stdlib.h>
#include "debug.h"
#include <math.h>

void dcheck_reproducibility_(
	int n, double* x, int incx, double* y, int incy,
	int* status, double* ref,
	I_double sum,
	I_double asum,
	I_double dot,
	I_double ref_nrm2,
	int order, int maxP) {
	// GENERATE DATA
	int i, j;
	double res;
	int P = 1;

	Idouble sum2;
	Idouble asum2;
	Idouble dot2;
	Idouble nrm2;
	double scale;
	Idouble tmp;
	int LN = (n + P - 1) / P;
	int lN;

	// PERMUTE DATA
	if (order == 0)
		dreverse(n, x);           // REVERSING
	else
		dsort_merge(n, x, order); // SORTING

	P = 1;
	while (P < n && P <= maxP) {
		if (P == 1)
			res = rdsum(n, x, incx);
		else {
			LN =  (n + P - 1) / P;
			dISetZero(sum2);
			for (j = 0; j < n; j += LN) {
				lN = LN < n - j ? LN : (n-j);
				tmp = dsumI(lN, x + j * incx, incx);
				// FLAT REDUCTION TREE
				dIAdd(&sum2, tmp);
			}
			res = Iconv2d(sum2);
		}
		if (res != ref[0]) {
			status[0]++;
			printf("\n SUM [P=%d,LN=%d]: %g %g %g ", P, LN, res, ref[0], res - ref[0]);
			if (P == 1) {
				sum2 = dsumI(n, x, incx);
			}
			printf("\n  :: ");
			dIprint(sum);
			printf("\n   : ");
			dIprint(sum2);
			return;
		}
		P *= 2;
	}

//	res = rdasum(n, x, 1);
	P = 1;
	while (P < n && P <= maxP) {
		if (P == 1)
			res = rdasum(n, x, incx);
		else {
			LN =  (n + P - 1) / P;
			dISetZero(asum2);
			for (j = 0; j < n; j += LN) {
				lN = LN < n - j ? LN : (n-j);
				tmp = dasumI(lN, x + j * incx, incx);
				// FLAT REDUCTION TREE
				dIAdd(&asum2, tmp);
			}
			res = Iconv2d(asum2);
		}
		P *= 2;
		if (res != ref[1]) {
			status[1]++;
			printf("\n ASUM : %g %g %g ", res, ref[1], res - ref[1]);

			if (P == 1) {
				asum2 = dasumI(n, x, incx);
			}
			printf("\n  :: ");
			dIprint(asum);
			printf("\n   : ");
			dIprint(asum2);
			return;
		}
	}

	P = 1;
	while (P < n && P <= maxP) {
		if (P == 1)
			res = rddot(n, x, incx, y, incy);
		else {
			LN =  (n + P - 1) / P;
			dISetZero(dot2);
			for (j = 0; j < n; j += LN) {
				lN = LN < n - j ? LN : (n-j);
				tmp = ddotI(lN, x + j * incx, incx, y + j * incy, incy);
				// FLAT REDUCTION TREE
				dIAdd(&dot2, tmp);
			}
			res = Iconv2d(dot2);
		}
		P *= 2;
		if (res != ref[2]) {
			status[2]++;
			printf("\n DOT : %g %g %g ", res, ref[2], res - ref[2]);

			if (P == 1) {
				dot2 = ddotI(n, x, incx, y, incy);
			}
			printf("\n  :: ");
			dIprint( dot);
			printf("\n   : ");
			dIprint( dot2);
			return;
		}
	}

	res = rdnrm2(n, x, incx);
	if (res != ref[3]) {
		status[3]++;
		printf("\n NRM2: %g %g %g ", res, ref[3], res - ref[3]);

		if (P == 1) {
			dISetZero(nrm2);
			scale = dnrm2I(n, x, incx, &nrm2);
		}
		printf("\n  :: ");
		dIprint( ref_nrm2);
		printf("\n   : ");
		dIprint( nrm2);
	}
}

void dcheck_reproducibility(int n, double* x, double* y,
	int* status, double* ref) {
	// GENERATE DATA
	int i, j;
	double res;
	int P = 1;
	for (i = 0;  i < 4; i++) status[i] = 0;

	Idouble sum;
	Idouble asum;
	Idouble dot;
	Idouble ref_nrm2;
	double scale;
	int maxP = 1024;

	// COMPUTE WITH ORIGINAL DATA

	sum  = dsumI(n, x, 1);
	asum = dasumI(n, x, 1);
	dot  = ddotI(n, x, 1, y, 1);

	dISetZero(ref_nrm2);
	scale = dnrm2I(n, x, 1, &ref_nrm2);

	ref[0] = Iconv2d(sum);
	ref[1] = Iconv2d(asum);
	ref[2] = Iconv2d(dot);
	ref[3] = scale * sqrt(Iconv2d(ref_nrm2));

	// REVERSE
	dcheck_reproducibility_(
		n, x, 1, y, 1, status, ref,
		sum, asum, dot, ref_nrm2,
		0, maxP);

	if ((status[0] + status[1] + status[2] + status[3]) > 0)
		return;

	for (i = 1; i < 3; i++) {
		dcheck_reproducibility_(
		n, x, 1, y, 1, status, ref,
		sum, asum, dot, ref_nrm2,
		i, maxP);
		if ((status[0] + status[1] + status[2] + status[3]) > 0)
			return;
		dcheck_reproducibility_(
		n, x, 1, y, 1, status, ref,
		sum, asum, dot, ref_nrm2,
		0, maxP);
		if (status[0] + status[1] + status[2] + status[3] > 0)
			return;
	}
	return;
}

