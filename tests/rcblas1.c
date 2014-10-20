#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <fenv.h>

#include <rblas.h>

#include "benchmark_macro.h"

int main( int argc, char **argv ) {
	float complex* v;
	float complex* y;
	float complex sum;
	int n;
	float ref;

	CHECK_HEADER_()

	// ALLOCATE MEMORY
	ref = 0.0;
	v = (float complex*)malloc(2*stop * incy * sizeof(float complex));
	y = (float complex*)malloc(2*stop * incy * sizeof(float complex));
	int i;
	for (i = 0; i < 2*stop*incy; i++) {
		CSET_(y[i], 1.0, 1.0);
	}

	if (print_header) {
	printf("\n");
	printf("%8s", "N");
	#ifdef CALL_SCASUM
	if (tests & (1 << ASUM_BIT))
	printf("%10s", "SCASUM");
	#endif
	if (tests & (1 << RASUM_BIT))
	printf("%10s", "RSCASUM");
	if (tests & (1 << RSUM_BIT))
	printf("%10s", "RZSUM");
	#ifdef CALL_SCNRM2
	if (tests & (1 << NRM2_BIT))
	printf("%10s", "SCNRM2");
	#endif
	if (tests & (1 << RNRM2_BIT))
	printf("%10s", "RSCNRM2");
	#ifdef CALL_CDOTC
	if (tests & (1 << DOT_BIT))
	printf("%10s", "CDOTC");
	#endif
	#ifdef CALL_CDOTU
	if (tests & (1 << DOT_BIT))
	printf("%10s", "CDOTU");
	#endif
	if (tests & (1 << RDOT_BIT))
	printf("%10s%10s", "RCDOTC","RCDOTU");

	printf("\n");
	}

	for (n = start; n < stop + step; n += step) {
		sgenvec_(n, (float*)v, 2*incv, dtype, 1.0);
		printf("%8d", n);

		//==== SCASUM ====
#		ifdef CALL_SCASUM
        if(flops){
          time1_(ref, tests, ASUM_BIT, CALL_SCASUM, n, v, incv, n * 4)
        }else{
          time1_(ref, tests, ASUM_BIT, CALL_SCASUM, n, v, incv, n)
        }
#		endif
		//==== RSCASUM ====
        if(flops){
          time_(ref, tests, RASUM_BIT, rscasum, n, v, incv, n * 22)
        }else{
          time_(ref, tests, RASUM_BIT, rscasum, n, v, incv, n)
        }
      

		//==== RDSUM ====
        if(flops){
          time_(sum, tests, RSUM_BIT, rcsum, n, v, incv, n * 20)
        }else{
          time_(sum, tests, RSUM_BIT, rcsum, n, v, incv, n)
        }

		//==== SCNRM2 ====
#		ifdef CALL_SCNRM2
        if(flops){
          time1_(ref, tests, NRM2_BIT, CALL_SCNRM2, n, v, incv, n * 10)
        }else{
          time1_(ref, tests, NRM2_BIT, CALL_SCNRM2, n, v, incv, n)
        }
#		endif

		//==== RSCNRM2 ====
        if(flops){
          time_(ref, tests, RNRM2_BIT, rscnrm2, n, v, incv, n * 24)
        }else{
          time_(ref, tests, RNRM2_BIT, rscnrm2, n, v, incv, n)
        }

		//==== DDOT ====
#		ifdef CALL_CDOTC
        if(flops){
          time3_(sum, tests, DOT_BIT, CALL_CDOTC, n, v, incv, y, incy, n * 6)
        }else{
          time3_(sum, tests, DOT_BIT, CALL_CDOTC, n, v, incv, y, incy, 2*n)
        }
#		endif
#		ifdef CALL_CDOTU
        if(flops){
          time3_(sum, tests, DOT_BIT, CALL_CDOTU, n, v, incv, y, incy, n * 6)
        }else{
          time3_(sum, tests, DOT_BIT, CALL_CDOTU, n, v, incv, y, incy, 2*n)
        }
#		endif

		//==== RDDOT ====
        if(flops){
          time2_(sum, tests, RDOT_BIT, rcdotc, n, v, incv, y, incy, 44*n)
        }else{
          time2_(sum, tests, RDOT_BIT, rcdotc, n, v, incv, y, incy, 2*n)
        }
        if(flops){
          time2_(sum, tests, RDOT_BIT, rcdotu, n, v, incv, y, incy, 44*n)
        }else{
          time2_(sum, tests, RDOT_BIT, rcdotu, n, v, incv, y, incy, 2*n)
        }
		fprintf(stdout, "\n");

	}

	free(v);
	free(y);
	return 0;
}
