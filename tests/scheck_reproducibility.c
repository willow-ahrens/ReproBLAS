#include <rblas.h>
#include <IndexedFP.h>
#include <stdio.h>
#include <stdlib.h>
#include "debug.h"
#include <math.h>

void scheck_reproducibility_(
	int n, float* x, int incx, float* y, int incy,
	int* status, float* ref,
	I_float sum,
	I_float asum,
	I_float dot,
	I_float ref_nrm2,
	int order, int maxP) {
	// GENERATE DATA
	int i, j;
	float res;
	int P = 1;

	Ifloat sum2;
	Ifloat asum2;
	Ifloat dot2;
	Ifloat nrm2;
	float scale;
	Ifloat tmp;
	int LN;
	int lN;

	// PERMUTE DATA
	if (order == 0)
		sreverse(n, x);           // REVERSING
	else
		ssort_merge(n, x, order); // SORTING

	P = 1;
	while (P < n && P <= maxP) {
		if (P == 1)
			res = rssum(n, x, 1);
		else {
			LN =  (n + P - 1) / P;
			sISetZero(sum2);
			for (j = 0; j < n; j += LN) {
				lN = LN < n - j ? LN : (n-j);
				tmp = ssumI(lN, x + j, 1);
				// FLAT REDUCTION TREE
//				printf("\nlocal sum: "); sIprint( tmp);
//				printf("\nsum: "); sIprint( sum2);
				sIAdd(&sum2, tmp);
//				printf("\nAggregation: "); sIprint( sum2);
			}
			res = Iconv2f(sum2);
		}
		if (res != ref[0]) {
			status[0]++;
			printf("\n SUM [P=%d]: %g %g %g ", P, res, ref[0], res - ref[0]);
			if (P == 1) {
				sum2 =  ssumI(n, x, 1);
			}
			printf("\n  :: ");
			sIprint( sum);
			printf("\n   : ");
			sIprint( sum2);
			return;
		}
		P *= 2;
	}

//	res = rsasum(n, x, 1);
	P = 1;
	while (P < n && P <= maxP) {
		if (P == 1)
			res = rsasum(n, x, 1);
		else {
			LN =  (n + P - 1) / P;
			sISetZero(asum2);
			for (j = 0; j < n; j += LN) {
				lN = LN < n - j ? LN : (n-j);
				tmp = sasumI(lN, x + j, 1);
				// FLAT REDUCTION TREE
				sIAdd(&asum2, tmp);
			}
			res = Iconv2f(asum2);
		}
		P *= 2;
		if (res != ref[1]) {
			status[1]++;
			printf("\n ASUM : %g %g %g ", res, ref[1], res - ref[1]);

			if (P == 1) {
				asum2 = sasumI(n, x, 1);
			}
			printf("\n  :: ");
			sIprint( asum);
			printf("\n   : ");
			sIprint( asum2);
			return;
		}
	}

	P = 1;
	while (P < n && P <= maxP) {
		if (P == 1)
			res = rsdot(n, x, 1, y, 1);
		else {
			LN =  (n + P - 1) / P;
			sISetZero(dot2);
			for (j = 0; j < n; j += LN) {
				lN = LN < n - j ? LN : (n-j);
				tmp = sdotI(lN, x + j, 1, y + j, 1);
				// FLAT REDUCTION TREE
				sIAdd(&dot2, tmp);
			}
			res = Iconv2f(dot2);
		}
		P *= 2;
		if (res != ref[2]) {
			status[2]++;
			printf("\n DOT : %g %g %g ", res, ref[2], res - ref[2]);

			if (P == 1) {
				dot2 = ssumI(n, x, 1);
			}
			printf("\n  :: ");
			sIprint( dot);
			printf("\n   : ");
			sIprint( dot2);
			return;
		}
	}

	res = rsnrm2(n, x, incx);
	if (res != ref[3]) {
		status[3]++;
		printf("\n NRM2: %g %g %g ", res, ref[3], res - ref[3]);

		if (P == 1) {
			sISetZero(nrm2);
			scale = snrm2I(n, x, 1, &nrm2);
		}
		printf("\n  :: ");
		sIprint( ref_nrm2);
		printf("\n   : ");
		sIprint( nrm2);
	}
}

void scheck_reproducibility(int n, float* x, float* y,
	int* status, float* ref) {
	// GENERATE DATA
	int i, j;
	float res;
	int P = 1;
	for (i = 0;  i < 4; i++) status[i] = 0;

	Ifloat sum;
	Ifloat asum;
	Ifloat dot;
	Ifloat ref_nrm2;
	float scale;
	int maxP = 1024;

	// COMPUTE WITH ORIGINAL DATA

	sum  = ssumI(n, x, 1);
	asum = sasumI(n, x, 1);
	dot  = sdotI(n, x, 1, y, 1);

	sISetZero(ref_nrm2);
	scale = snrm2I(n, x, 1, &ref_nrm2);

	ref[0] = Iconv2f(sum);
	ref[1] = Iconv2f(asum);
	ref[2] = Iconv2f(dot);
	ref[3] = scale * sqrt(Iconv2f(ref_nrm2));

	// REVERSE
	scheck_reproducibility_(
		n, x, 1, y, 1, status, ref,
		sum, asum, dot, ref_nrm2,
		0, maxP);

	if ((status[0] + status[1] + status[2] + status[3]) > 0)
		return;

	for (i = 1; i < 3; i++) {
		scheck_reproducibility_(
		n, x, 1, y, 1, status, ref,
		sum, asum, dot, ref_nrm2,
		i, maxP);
		if ((status[0] + status[1] + status[2] + status[3]) > 0)
			return;
		scheck_reproducibility_(
		n, x, 1, y, 1, status, ref,
		sum, asum, dot, ref_nrm2,
		0, maxP);
		if ((status[0] + status[1] + status[2] + status[3]) > 0)
			return;
	}
	return;
}


