#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <fenv.h>

#include <rblas.h>

#include "benchmark_macro.h"

int main( int argc, char **argv ) {
	double complex* v;
	double complex* y;
	double complex sum;
	int n;
	double ref;

	CHECK_HEADER_()

	// ALLOCATE MEMORY
	ref = 0.0;
	v = (double complex*)malloc(2*stop * sizeof(double complex));
	y = (double complex*)malloc(2*stop * sizeof(double complex));
	int i;
	for (i = 0; i < 2*stop; i++) {
		ZSET_(y[i], 1.0, 1.0);
	}

	if (print_header) {
	printf("\n");
	printf("%8s", "N");
	#ifdef CALL_DZASUM
	if (tests & (1 << ASUM_BIT))
	printf("%10s", "DZASUM");
	#endif
	if (tests & (1 << RASUM_BIT))
	printf("%10s", "RDZASUM");
	if (tests & (1 << RSUM_BIT))
	printf("%10s", "RZSUM");
	#ifdef CALL_DZNRM2
	if (tests & (1 << NRM2_BIT))
	printf("%10s", "DZNRM2");
	#endif
	if (tests & (1 << RNRM2_BIT))
	printf("%10s", "RDZNRM2");
	#ifdef CALL_ZDOTC
	if (tests & (1 << DOT_BIT))
	printf("%10s", "ZDOTC");
	#endif
	#ifdef CALL_ZDOTU
	if (tests & (1 << DOT_BIT))
	printf("%10s", "ZDOTU");
	#endif
	if (tests & (1 << RDOT_BIT))
	printf("%10s%10s", "RZDOTC","RZDOTU");

	printf("\n");
	}

	for (n = start; n < stop + step; n += step) {
		dgenvec_(n, (double*)v, 2*incv, dtype, 1.0);
		printf("%8d", n);

		//==== DZASUM ====
#		ifdef CALL_DZASUM
		time1_(ref, tests, ASUM_BIT, CALL_DZASUM, n, v, incv, n)
#		endif
		//==== RDZASUM ====
		time_(ref,tests, RASUM_BIT, rdzasum, n, v, incv, n)

		//==== RDSUM ====
		time_(sum,tests, RSUM_BIT, rzsum, n, v, incv, n)

		//==== DZNRM2 ====
#		ifdef CALL_DZNRM2
		time1_(ref, tests, NRM2_BIT, CALL_DZNRM2, n, v, incv, n)
#		endif

		//==== RDZNRM2 ====
		time_(ref,tests, RNRM2_BIT, rdznrm2, n, v, incv, n)

		//==== DDOT ====
#		ifdef CALL_ZDOTC
		time3_(sum, tests, DOT_BIT, CALL_ZDOTC, n, v, incv, y, incy, 2*n)
#		endif
#		ifdef CALL_ZDOTU
		time3_(sum, tests, DOT_BIT, CALL_ZDOTU, n, v, incv, y, incy, 2*n)
#		endif

		//==== RDDOT ====
		time2_(sum, tests, RDOT_BIT, rzdotc, n, v, incv, y, incy, 2*n)
		time2_(sum, tests, RDOT_BIT, rzdotu, n, v, incv, y, incy, 2*n)

		fprintf(stdout, "\n");

	}

	free(v);
	free(y);
	return 0;
}
