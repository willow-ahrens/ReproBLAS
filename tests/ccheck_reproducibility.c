#include <rblas.h>
#include <stdio.h>
#include <stdlib.h>
#include "debug.h"
#include <math.h>

#define zeq(X,Y) ((CREAL_(X) == CREAL_(Y)) && (CIMAG_(X) == CIMAG_(Y)))
#define cneq(X,Y) ((CREAL_(X) != CREAL_(Y)) || (CIMAG_(X) != CIMAG_(Y)))

void ccheck_reproducibility_(
	int n, float complex* x, int incx, float complex* y, int incy,
	int* status, float complex* ref,
	I_float_Complex sum,
	I_float asum,
	I_float_Complex dotu,
	I_float_Complex dotc,
	I_float ref_nrm2,
	int order, int maxP) {
	// GENERATE DATA
	int i, j;
	float complex res;
	float dres;
	int P = 1;

	I_float_Complex sum2;
	Ifloat asum2;
	I_float_Complex dotu2;
	I_float_Complex dotc2;
	Ifloat nrm2;
	float scale;
	Ifloat tmp;
	I_float_Complex ztmp;
	int LN = (n + P - 1) / P;
	int lN;

	// PERMUTE DATA
	if (order == 0)
		creverse(n, x, incx);           // REVERSING
	else
		csort_merge(n, x, incx, order); // SORTING

	P = 1;
	while (P < n && P <= maxP) {
		if (P == 1)
			res = rcsum(n, x, incx);
		else {
			LN =  (n + P - 1) / P;
			cISetZero(sum2);
			for (j = 0; j < n; j += LN) {
				lN = LN < n - j ? LN : (n-j);
				ztmp = csumI(lN, x + j * incx, incx);
				// FLAT REDUCTION TREE
				cIAdd(&sum2, ztmp);
			}
			res = Iconv2c(sum2);
		}
		if (cneq(res, ref[0])) {
			status[0]++;
			printf("\n SUM [P=%d,LN=%d]: {%g , %g} - {%g , %g} = {%g , %g} ", P, LN,
				ZREAL_(res), ZIMAG_(res), ZREAL_(ref[0]), ZIMAG_(ref[0]),
				ZREAL_(res) - ZREAL_(ref[0]), ZIMAG_(ref[0]) - ZIMAG_(res));
			if (P == 1) {
				sum2 = csumI(n, x, incx);
			}
			printf("\n  :: ");
			cIprint(sum);
			printf("\n   : ");
			cIprint(sum2);
			return;
		}
		P *= 2;
	}

//	res = rdasum(n, x, 1);
	P = 1;
	while (P < n && P <= maxP) {
		if (P == 1)
			dres = rscasum(n, x, incx);
		else {
			LN =  (n + P - 1) / P;
			dISetZero(asum2);
			for (j = 0; j < n; j += LN) {
				lN = LN < n - j ? LN : (n-j);
				tmp = scasumI(lN, x + j*incx, incx);
				// FLAT REDUCTION TREE
				sIAdd(&asum2, tmp);
			}
			dres = Iconv2f(asum2);
		}
		CSET_(res, dres, 0);
		if (cneq(res, ref[1])) {
			status[1]++;
//			printf("\n ASUM : %g %g %g ", res, ref[1], res - ref[1]);
			printf("\nASUM [P=%d,LN=%d]: \n", P, LN);

			if (P == 1) {
				asum2 = scasumI(n, x, incx);
			}
			printf("\n  :: ");
			sIprint(asum);
			printf("\n   : ");
			sIprint(asum2);
			return;
		}
		P *= 2;
	}

	P = 1;
	while (P < n && P <= maxP) {
		if (P == 1)
			res = rcdotu(n, x, incx, y, incy);
		else {
			LN =  (n + P - 1) / P;
			cISetZero(dotu2);
			for (j = 0; j < n; j += LN) {
				lN = LN < n - j ? LN : (n-j);
				ztmp = cdotuI(lN, x + j*incx, incx, y + j*incy, incy);
				// FLAT REDUCTION TREE
				cIAdd(&dotu2, ztmp);
			}
			res = Iconv2c(dotu2);
		}
		if (cneq(res, ref[2])) {
			status[2]++;
			printf("\n DOTU [P=%d,LN=%d]: {%g , %g} - {%g , %g} = {%g , %g} ", P, LN,
				ZREAL_(res), ZIMAG_(res), ZREAL_(ref[2]), ZIMAG_(ref[2]),
				ZREAL_(res) - ZREAL_(ref[2]), ZIMAG_(ref[2]) - ZIMAG_(res));

			if (P == 1) {
				dotu2 = cdotuI(n, x, incx, y, incy);
			}
			printf("\n  :: ");
			cIprint( dotu );
			printf("\n   : ");
			cIprint( dotu2);
			return;
		}
		P *= 2;
	}

	P = 1;
	while (P < n && P <= maxP) {
		if (P == 1)
			res = rcdotc(n, x, incx, y, incy);
		else {
			LN =  (n + P - 1) / P;
			cISetZero(dotc2);
			for (j = 0; j < n; j += LN) {
				lN = LN < n - j ? LN : (n-j);
				ztmp = cdotcI(lN, x + j*incx, incx, y + j*incy, incy);
				// FLAT REDUCTION TREE
				cIAdd(&dotc2, ztmp);
			}
			res = Iconv2c(dotc2);
		}
		if (cneq(res, ref[3])) {
			status[3]++;
			printf("\n DOTC [P=%d,LN=%d]: {%g , %g} - {%g , %g} = {%g , %g} ", P, LN,
				ZREAL_(res), ZIMAG_(res), ZREAL_(ref[3]), ZIMAG_(ref[3]),
				ZREAL_(res) - ZREAL_(ref[3]), ZIMAG_(ref[3]) - ZIMAG_(res));

			if (P == 1) {
				dotu2 = cdotcI(n, x, incx, y, incy);
			}
			printf("\n  :: ");
			cIprint( dotc );
			printf("\n   : ");
			cIprint( dotc2);
			return;
		}
		P *= 2;
	}

	CSET_(res, rscnrm2(n, x, incx), 0);
	if (cneq(res, ref[4])) {
		status[4]++;
//		printf("\n NRM2: %g %g %g ", res, ref[3], res - ref[3]);

		if (P == 1) {
			dISetZero(nrm2);
			scale = scnrm2I(n, x, incx, &nrm2);
		}
		printf("\n  :: ");
		sIprint( ref_nrm2);
		printf("\n   : ");
		sIprint( nrm2);
	}
}

void ccheck_reproducibility(int n, float complex* x, int incx, float complex* y, int incy,
	int* status, float complex* ref) {
	// GENERATE DATA
	int i, j;
	float res;
	int P = 1;
	for (i = 0;  i < 5; i++) status[i] = 0;

	I_float_Complex sum;
	Ifloat asum;
	I_float_Complex dotu;
	I_float_Complex dotc;
	Ifloat ref_nrm2;
	float scale;
	int maxP = 1024;

	// COMPUTE WITH ORIGINAL DATA

	sum  = csumI(n, x, incx);
	asum = scasumI(n, x, incx);
	dotu = cdotuI(n, x, incx, y, incy);
	dotc = cdotcI(n, x, incx, y, incy);

	dISetZero(ref_nrm2);
	scale = scnrm2I(n, x, incx, &ref_nrm2);

	ref[0] = Iconv2c(sum);
	CSET_(ref[1], Iconv2f(asum), 0);
	ref[2] = Iconv2c(dotu);
	ref[3] = Iconv2c(dotc);
	CSET_(ref[4], scale * sqrt(Iconv2f(ref_nrm2)), 0);

	// REVERSE
	ccheck_reproducibility_(
		n, x, incx, y, incy, status, ref,
		sum, asum, dotu, dotc, ref_nrm2,
		0, maxP);

	if ((status[0] + status[1] + status[2] + status[3]) > 0) {
		printf("FAILED AT REVERSING \n");
		return;
	}

	for (i = 1; i < 3; i++) {
		ccheck_reproducibility_(
		n, x, incx, y, incy, status, ref,
		sum, asum, dotu, dotc, ref_nrm2,
		i, maxP);
		if ((status[0] + status[1] + status[2] + status[3] + status[4]) > 0) {
			printf("FAILED AT ORDER: %d\n", i);
			return;
		}
		ccheck_reproducibility_(
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

/*
void ccheck_reproducibility(int n, float complex* x, int incx, float complex* y, int incy,
	int* status, float complex* ref) {
	// GENERATE DATA
	int i;
	float complex res;

	// original data
	Ifloat   ref_asum;
	Ifloat   asum;
	sIcomplex ref_sum;
	sIcomplex ref_dotc;
	sIcomplex ref_dotu;
	sIcomplex sum;
	sIcomplex dotc;
	sIcomplex dotu;
	float   ref_nrm2[10];
	float   nrm2[10];

	sISetZero(ref_asum);
	sISetZero(ref_nrm2 + 1);
	cISetZero(ref_sum);
	cISetZero(ref_dotc);
	cISetZero(ref_dotu);

	csumI  (n, x, incx, ref_sum);
	scasumI(n, x, incx, ref_asum);
	cdotcI (n, x, incx, y, incy, ref_dotc);
	cdotuI (n, x, incx, y, incy, ref_dotu);
	ref_nrm2[0] = scnrm2I(n, x, incx, ref_nrm2 + 1);

	ref[0] = Iconv2c(ref_sum);
	CSET_(ref[1], Iconv2f(ref_asum), 0.0);
	ref[2]= Iconv2c(ref_dotc);
	ref[3]= Iconv2c(ref_dotu);
	CSET_(ref[4], ref_nrm2[0] * sqrt(Iconv2f(ref_nrm2 + 1)), 0.0);

	// reverse
	creverse(n, x, incx);

	res = rcsum(n, x, incx);
	if cneq(res , ref[0]) {
		status[0]++;
		printf("\n SUM reverse: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[0]),
			CIMAG_(res) - CIMAG_(ref[0]));
//		for (i = 0; i < 6; i++) sum2[i] = 0.0;
//		dsumI(n, x, 1, 3, sum2);
//		printf("\n  :: ");
//		sIprint(3, sum);
//		printf("\n   : ");
//		sIprint(3, sum2);
	}

	CSET_(res, rscasum(n, x, incx), 0.0);
	if cneq(res , ref[1]) {
		status[1]++;
		printf("\n ASUM reverse: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[1]),
			CIMAG_(res) - CIMAG_(ref[1]));
//		for (i = 0; i < 6; i++) asum2[i] = 0.0;
//		dsumI(n, x, 1, 3, asum2);
//		printf("\n  :: ");
//		sIprint(3, asum);
//		printf("\n   : ");
//		sIprint(3, asum2);
	}

	res = rcdotc(n, x, incx, y, incy);
	if cneq(res , ref[2]) {
		status[2]++;
		printf("\n DOTC reverse: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[2]),
			CIMAG_(res) - CIMAG_(ref[2]));
	}


	res = rcdotu(n, x, incx, y, incy);
	if cneq(res , ref[3]) {
		status[3]++;
		printf("\n DOTU reverse: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[3]),
			CIMAG_(res) - CIMAG_(ref[3]));
	}

	CSET_(res, rscnrm2(n, x, incx), 0.0);
//	printf("nrm2 = %g \n", CREAL_(res));
	if cneq(res , ref[4]) {
		status[4]++;
		printf("\n NRM2 reverse: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[4]),
			CIMAG_(res) - CIMAG_(ref[4]));
//		for (i = 0; i < 6; i++) nrm2[i] = 0.0;
//		nrm2[6] = dnrm2I(n, x, 1, 3, nrm2);
//		printf("\n  :: ");
//		sIprint(3, ref_nrm2);
//		printf("\n   : ");
//		sIprint(3, nrm2);
	}

	//=====================================
	// INCREASING ORDER
	csort_merge(n, x, incx, 1);
//	csort_bubble(n, x, incx, 1);

	res = rcsum(n, x, incx);
	if cneq(res , ref[0]) {
		status[0]++;
		printf("\n SUM +++: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[0]),
			CIMAG_(res) - CIMAG_(ref[0]));
//		for (i = 0; i < 6; i++) sum2[i] = 0.0;
//		dsumI(n, x, 1, 3, sum2);
//		printf("\n  :: ");
//		sIprint(3, sum);
//		printf("\n   : ");
//		sIprint(3, sum2);
	}

	CSET_(res, rscasum(n, x, incx),0.0);
	if cneq(res , ref[1]) {
		status[1]++;
		printf("\n ASUM +++: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[1]),
			CIMAG_(res) - CIMAG_(ref[1]));
//		for (i = 0; i < 6; i++) asum2[i] = 0.0;
//		dsumI(n, x, 1, 3, asum2);
//		printf("\n  :: ");
//		sIprint(3, asum);
//		printf("\n   : ");
//		sIprint(3, asum2);
	}

	res = rcdotc(n, x, incx, y, incy);
	if cneq(res , ref[2]) {
		status[2]++;
		printf("\n DOTC +++: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[2]),
			CIMAG_(res) - CIMAG_(ref[2]));
	}

	res = rcdotu(n, x, incx, y, incy);
	if cneq(res , ref[3]) {
		status[3]++;
		printf("\n DOTU +++: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[3]),
			CIMAG_(res) - CIMAG_(ref[3]));
	}
	CSET_(res, rscnrm2(n, x, incx), 0.0);
	if cneq(res , ref[4]) {
		status[4]++;
		printf("\n NRM2 +++: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[4]),
			CIMAG_(res) - CIMAG_(ref[4]));
//		for (i = 0; i < 6; i++) nrm2[i] = 0.0;
//		dnrm2I(n, x, 1, 3, nrm2);
//		printf("\n  :: ");
//		sIprint(3, ref_nrm2);
//		printf("\n   : ");
//		sIprint(3, nrm2);
	}

	//=====================================
	// DECREASING ORDER
	creverse(n, x, incx);

	res = rcsum(n, x, incx);
	if cneq(res , ref[0]) {
		status[0]++;
		printf("\n SUM ---: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[0]),
			CIMAG_(res) - CIMAG_(ref[0]));
//		for (i = 0; i < 6; i++) sum2[i] = 0.0;
//		dsumI(n, x, 1, 3, sum2);
//		printf("\n  :: ");
//		sIprint(3, sum);
//		printf("\n   : ");
//		sIprint(3, sum2);
	}

	CSET_(res, rscasum(n, x, incx), 0.0);
	if cneq(res , ref[1]) {
		status[1]++;
		printf("\n ASUM ---: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[1]),
			CIMAG_(res) - CIMAG_(ref[1]));
//		for (i = 0; i < 6; i++) asum2[i] = 0.0;
//		dsumI(n, x, 1, 3, asum2);
//		printf("\n  :: ");
//		sIprint(3, asum);
//		printf("\n   : ");
//		sIprint(3, asum2);
	}

	res = rcdotc(n, x, incx, y, incy);
	if cneq(res , ref[2]) {
		status[2]++;
		printf("\n DOTC ---: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[2]),
			CIMAG_(res) - CIMAG_(ref[2]));
	}

	res = rcdotu(n, x, incx, y, incy);
	if cneq(res , ref[3]) {
		status[3]++;
		printf("\n DOTU ---: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[3]),
			CIMAG_(res) - CIMAG_(ref[3]));
	}

	CSET_(res, rscnrm2(n, x, incx), 0.0);
	if cneq(res , ref[4]) {
		status[4]++;
		printf("\n NRM2 ---: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[4]),
			CIMAG_(res) - CIMAG_(ref[4]));
//		for (i = 0; i < 6; i++) nrm2[i] = 0.0;
//		dnrm2I(n, x, 1, 3, nrm2);
//		printf("\n  :: ");
//		sIprint(3, ref_nrm2);
//		printf("\n   : ");
//		sIprint(3, nrm2);
	}

	//=====================================
	// INCREASING MAGNITUDE
	csort_merge(n, x, incx, 2);
//	csort_bubble(n, x, incx, 2);
	res = rcsum(n, x, incx);
	if cneq(res , ref[0]) {
		status[0]++;
		printf("\n SUM ++MAG: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[0]),
			CIMAG_(res) - CIMAG_(ref[0]));
//		for (i = 0; i < 6; i++) sum2[i] = 0.0;
//		dsumI(n, x, 1, 3, sum2);
//		printf("\n  :: ");
//		sIprint(3, sum);
//		printf("\n   : ");
//		sIprint(3, sum2);
	}

	CSET_(res, rscasum(n, x, incx), 0.0);
	if cneq(res , ref[1]) {
		status[1]++;
		printf("\n ASUM ++MAG: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[1]),
			CIMAG_(res) - CIMAG_(ref[1]));
//		for (i = 0; i < 6; i++) asum2[i] = 0.0;
//		dsumI(n, x, 1, 3, asum2);
//		printf("\n  :: ");
//		sIprint(3, asum);
//		printf("\n   : ");
//		sIprint(3, asum2);
	}

	res = rcdotc(n, x, incx, y, incy);
	if cneq(res , ref[2]) {
		status[2]++;
		printf("\n DOTC ++MAG: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[2]),
			CIMAG_(res) - CIMAG_(ref[2]));
	}

	res = rcdotu(n, x, incx, y, incy);
	if cneq(res , ref[3]) {
		status[3]++;
		printf("\n DOTU ++MAG: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[3]),
			CIMAG_(res) - CIMAG_(ref[3]));
	}

	CSET_(res, rscnrm2(n, x, incx), 0.0);
	if cneq(res , ref[4]) {
		status[4]++;
		printf("\n NRM2 ++MAG: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[4]),
			CIMAG_(res) - CIMAG_(ref[4]));
//		for (i = 0; i < 6; i++) nrm2[i] = 0.0;
//		dnrm2I(n, x, 1, 3, nrm2);
//		printf("\n  :: ");
//		sIprint(3, ref_nrm2);
//		printf("\n   : ");
//		sIprint(3, nrm2);
	}

	//=====================================
	// DECREASING MAGNITUDE
	creverse(n, x, incx);

	res = rcsum(n, x, incx);
	if cneq(res , ref[0]) {
		status[0]++;
		printf("\n SUM --MAG: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[0]),
			CIMAG_(res) - CIMAG_(ref[0]));
//		for (i = 0; i < 6; i++) asum2[i] = 0.0;
//		dsumI(n, x, 1, 3, asum2);
//		printf("\n  :: ");
//		sIprint(3, asum);
//		printf("\n   : ");
//		sIprint(3, asum2);
	}

	CSET_(res, rscasum(n, x, incx), 0.0);
	if cneq(res , ref[1]) {
		status[1]++;
		printf("\n ASUM --MAG: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[1]),
			CIMAG_(res) - CIMAG_(ref[1]));
	}

	res = rcdotc(n, x, incx, y, incy);
	if cneq(res , ref[2]) {
		status[2]++;
		printf("\n DOTC --MAG: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[2]),
			CIMAG_(res) - CIMAG_(ref[2]));
	}

	res = rcdotu(n, x, incx, y, incy);
	if cneq(res , ref[3]) {
		status[3]++;
		printf("\n DOTU --MAG: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[3]),
			CIMAG_(res) - CIMAG_(ref[3]));
	}
	
	CSET_(res, rscnrm2(n, x, incx), 0.0);
	if cneq(res , ref[4]) {
		status[4]++;
		printf("\n NRM2 --MAG: %g %g %g %g ",
			CREAL_(res), CIMAG_(res),
			CREAL_(res) - CREAL_(ref[4]),
			CIMAG_(res) - CIMAG_(ref[4]));
//		for (i = 0; i < 6; i++) nrm2[i] = 0.0;
//		dnrm2I(n, x, 1, 3, nrm2);
//		printf("\n  :: ");
//		sIprint(3, ref_nrm2);
//		printf("\n   : ");
//		sIprint(3, nrm2);
	}

}
*/
