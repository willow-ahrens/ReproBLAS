#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <fenv.h>

#include <rblas.h>

#include "benchmark_macro.h"

int main( int argc, char **argv ) {
	float* v;
	float* y;
	float sum;
	float lsum[32];
	int n;
	float ref;

	CHECK_HEADER_()

	// ALLOCATE MEMORY
	ref = 0.0;
	v = (float*)malloc(stop * incv * sizeof(float));
	y = (float*)malloc(stop * incy * sizeof(float));
	
	int i;
	for (i = 0; i < stop*incy; i++)
		y[i] = 1.0;

	if (print_header) {
	printf("\n");
	printf("%8s", "N");
	#ifdef CALL_DASUM
	if (tests & (1 << ASUM_BIT))
	printf("%10s", "SASUM");
	#endif
	if (tests & (1 << RASUM_BIT))
	printf("%10s", "RSASUM");
	if (tests & (1 << RSUM_BIT))
	printf("%10s", "RSSUM");
	#ifdef CALL_DNRM2
	if (tests & (1 << NRM2_BIT))
	printf("%10s", "SNRM2");
	#endif
	if (tests & (1 << RNRM2_BIT))
	printf("%10s", "RSNRM2");
	#ifdef CALL_DDOT
	if (tests & (1 << DOT_BIT))
	printf("%10s", "SDOT");
	#endif
	if (tests & (1 << RDOT_BIT))
	printf("%10s", "RSDOT");

	printf("\n");
	}

	for (n = start; n < stop + step; n += step) {
		sgenvec_(n, v, incv, dtype, 1.0);
		printf("%8d", n);

		//==== DASUM ====
#		ifdef CALL_DASUM
		time1_(ref, tests, ASUM_BIT, CALL_SASUM, n, v, incv, n)
#		endif
		//==== RDASUM ====
		time_(sum,tests, RASUM_BIT, rsasum, n, v, incv, n)

		//==== RDSUM ====
		time_(sum,tests, RSUM_BIT, rssum, n, v, incv, n)

		//==== DNRM2 ====
#		ifdef CALL_DNRM2
		time1_(ref, tests, NRM2_BIT, CALL_SNRM2, n, v, incv, n)
#		endif

		//==== RDNRM2 ====
		time_(sum,tests, RNRM2_BIT, rsnrm2, n, v, incv, n)

		//==== DDOT ====
#		ifdef CALL_DDOT
		time3_(ref, tests, DOT_BIT, CALL_SDOT, n, v, incv, y, incy, 2*n)
#		endif

		//==== RDDOT ====
		time2_(sum, tests, RDOT_BIT, rsdot, n, v, incv, y, incy, 2*n)

		fprintf(stdout, "\n");

	}

	free(v);
	free(y);
	return 0;
}
