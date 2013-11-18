#include <rblas.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "debug.h"

#define zeq(X,Y) ((ZREAL(X) == ZREAL_(Y)) && (ZIMAG_(X) == ZIMAG_(Y)))
#define zneq(X,Y) ((ZREAL_(X) != ZREAL_(Y)) || (ZIMAG_(X) != ZIMAG_(Y)))

void zcheck_reproducibility_(
	int n, double complex* x, int incx, double complex* y, int incy,
	int* status, double complex* ref,
	I_double_Complex sum,
	I_double asum,
	I_double_Complex dotu,
	I_double_Complex dotc,
	I_double ref_nrm2,
	int order, int maxP) {
	// GENERATE DATA
	int i, j;
	double complex res;
	double dres;
	int P = 1;

	I_double_Complex sum2;
	Idouble asum2;
	I_double_Complex dotu2;
	I_double_Complex dotc2;
	Idouble nrm2;
	double scale;
	Idouble tmp;
	I_double_Complex ztmp;
	int LN = (n + P - 1) / P;
	int lN;

	// PERMUTE DATA
	if (order == 0)
		zreverse(n, x, incx);           // REVERSING
	else
		zsort_merge(n, x, incx, order); // SORTING

	P = 1;
	while (P < n && P <= maxP) {
		if (P == 1)
			res = rzsum(n, x, incx);
		else {
			LN =  (n + P - 1) / P;
			zISetZero(sum2);
			for (j = 0; j < n; j += LN) {
				lN = LN < n - j ? LN : (n-j);
				ztmp = zsumI(lN, x + j * incx, incx);
				// FLAT REDUCTION TREE
				zIAdd(&sum2, ztmp);
			}
			res = Iconv2z(sum2);
		}
		if (zneq(res, ref[0])) {
			status[0]++;
			printf("\n SUM [P=%d,LN=%d]: {%g , %g} - {%g , %g} = {%g , %g} ", P, LN,
				ZREAL_(res), ZIMAG_(res), ZREAL_(ref[0]), ZIMAG_(ref[0]),
				ZREAL_(res) - ZREAL_(ref[0]), ZIMAG_(ref[0]) - ZIMAG_(res));
			if (P == 1) {
				sum2 = zsumI(n, x, incx);
			}
			printf("\n  :: ");
			zIprint(sum);
			printf("\n   : ");
			zIprint(sum2);
			return;
		}
		P *= 2;
	}

//	res = rdasum(n, x, 1);
	P = 1;
	while (P < n && P <= maxP) {
		if (P == 1)
			dres = rdzasum(n, x, incx);
		else {
			LN =  (n + P - 1) / P;
			dISetZero(asum2);
			for (j = 0; j < n; j += LN) {
				lN = LN < n - j ? LN : (n-j);
				tmp = dzasumI(lN, x + j*incx, incx);
				// FLAT REDUCTION TREE
				dIAdd(&asum2, tmp);
			}
			dres = Iconv2d(asum2);
		}
		ZSET_(res, dres, 0);
		if (zneq(res, ref[1])) {
			status[1]++;
//			printf("\n ASUM : %g %g %g ", res, ref[1], res - ref[1]);
			printf("\nASUM [P=%d,LN=%d]: \n", P, LN);

			if (P == 1) {
				asum2 = dzasumI(n, x, incx);
			}
			printf("\n  :: ");
			dIprint(asum);
			printf("\n   : ");
			dIprint(asum2);
			return;
		}
		P *= 2;
	}

	P = 1;
	while (P < n && P <= maxP) {
		if (P == 1)
			res = rzdotu(n, x, incx, y, incy);
		else {
			LN =  (n + P - 1) / P;
			zISetZero(dotu2);
			for (j = 0; j < n; j += LN) {
				lN = LN < n - j ? LN : (n-j);
				ztmp = zdotuI(lN, x + j*incx, incx, y + j*incy, incy);
				// FLAT REDUCTION TREE
				zIAdd(&dotu2, ztmp);
			}
			res = Iconv2z(dotu2);
		}
		if (zneq(res, ref[2])) {
			status[2]++;
			printf("\n DOTU [P=%d,LN=%d]: {%g , %g} - {%g , %g} = {%g , %g} ", P, LN,
				ZREAL_(res), ZIMAG_(res), ZREAL_(ref[2]), ZIMAG_(ref[2]),
				ZREAL_(res) - ZREAL_(ref[2]), ZIMAG_(ref[2]) - ZIMAG_(res));

			if (P == 1) {
				dotu2 = zdotuI(n, x, incx, y, incy);
			}
			printf("\n  :: ");
			zIprint( dotu );
			printf("\n   : ");
			zIprint( dotu2);
			return;
		}
		P *= 2;
	}

	P = 1;
	while (P < n && P <= maxP) {
		if (P == 1)
			res = rzdotc(n, x, incx, y, incy);
		else {
			LN =  (n + P - 1) / P;
			zISetZero(dotc2);
			for (j = 0; j < n; j += LN) {
				lN = LN < n - j ? LN : (n-j);
				ztmp = zdotcI(lN, x + j*incx, incx, y + j*incy, incy);
				// FLAT REDUCTION TREE
				zIAdd(&dotc2, ztmp);
			}
			res = Iconv2z(dotc2);
		}
		if (zneq(res, ref[3])) {
			status[3]++;
			printf("\n DOTC [P=%d,LN=%d]: {%g , %g} - {%g , %g} = {%g , %g} ", P, LN,
				ZREAL_(res), ZIMAG_(res), ZREAL_(ref[3]), ZIMAG_(ref[3]),
				ZREAL_(res) - ZREAL_(ref[3]), ZIMAG_(ref[3]) - ZIMAG_(res));

			if (P == 1) {
				dotu2 = zdotcI(n, x, incx, y, incy);
			}
			printf("\n  :: ");
			zIprint( dotc );
			printf("\n   : ");
			zIprint( dotc2);
			return;
		}
		P *= 2;
	}

	ZSET_(res, rdznrm2(n, x, incx), 0);
	if (zneq(res, ref[4])) {
		status[4]++;
//		printf("\n NRM2: %g %g %g ", res, ref[3], res - ref[3]);

		if (P == 1) {
			dISetZero(nrm2);
			scale = dznrm2I(n, x, incx, &nrm2);
		}
		printf("\n  :: ");
		dIprint( ref_nrm2);
		printf("\n   : ");
		dIprint( nrm2);
	}
}

void zcheck_reproducibility(int n, double complex* x, int incx, double complex* y, int incy,
	int* status, double complex* ref) {
	// GENERATE DATA
	int i, j;
	double res;
	int P = 1;
	for (i = 0;  i < 5; i++) status[i] = 0;

	I_double_Complex sum;
	Idouble asum;
	I_double_Complex dotu;
	I_double_Complex dotc;
	Idouble ref_nrm2;
	double scale;
	int maxP = 1024;

	// COMPUTE WITH ORIGINAL DATA

	sum  = zsumI(n, x, incx);
	asum = dzasumI(n, x, incx);
	dotu = zdotuI(n, x, incx, y, incy);
	dotc = zdotcI(n, x, incx, y, incy);

	dISetZero(ref_nrm2);
	scale = dznrm2I(n, x, incx, &ref_nrm2);

	ref[0] = Iconv2z(sum);
	ZSET_(ref[1], Iconv2d(asum), 0);
	ref[2] = Iconv2z(dotu);
	ref[3] = Iconv2z(dotc);
	ZSET_(ref[4], scale * sqrt(Iconv2d(ref_nrm2)), 0);

	// REVERSE
	zcheck_reproducibility_(
		n, x, incx, y, incy, status, ref,
		sum, asum, dotu, dotc, ref_nrm2,
		0, maxP);

	if ((status[0] + status[1] + status[2] + status[3]) > 0) {
		printf("FAILED AT REVERSING \n");
		return;
	}

	for (i = 1; i < 3; i++) {
		zcheck_reproducibility_(
		n, x, incx, y, incy, status, ref,
		sum, asum, dotu, dotc, ref_nrm2,
		i, maxP);
		if ((status[0] + status[1] + status[2] + status[3] + status[4]) > 0) {
			printf("FAILED AT ORDER: %d\n", i);
			return;
		}
		zcheck_reproducibility_(
		n, x, incx, y, incy, status, ref,
		sum, asum, dotu, dotc, ref_nrm2,
		0, maxP);
		if (status[0] + status[1] + status[2] + status[3] + status[4] > 0) {
			printf("FAILED AT ORDER: %d\n", -i);
			return;
		}
	}
	return;
}

