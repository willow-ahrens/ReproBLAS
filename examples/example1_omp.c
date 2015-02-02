#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include <IndexedFP.h>
#include "doubledouble.h"

//double v(int i, int n) {
//	return  sin(2*M_PI * (i / (double)n - 0.5));
//}
#define v(i,n) sin(2 * M_PI * (i / (double)n - 0.5))

int main (int argc, char** args) {
	int n = 1000000;
	int i;
	double s1, s2;
	double ss[4];
	I_double I_s;
	int chunk = 128;

	if (argc > 1) {
		int nbthreads = atoi(args[1]);
		omp_set_num_threads(nbthreads);
	}

	printf("\nSum of sin(2 * M_PI * (i / n - 0.5).    n = %d \n",n);
	// SUMMATION USING PRIMITIVE TYPE
	s1 = 0;
	for (i = 0; i < n; i++) {
		s1 += v(i,n);
	}
	printf("\n>>%20s :  %.17g \n", "double [sequential]", s1);
	// multi-threaded
	s2 = 0;
	#pragma omp parallel
	{
		if (omp_get_thread_num() == 0)
			printf("           %2d theads ", omp_get_num_threads());
		#pragma omp for reduction(+:s2)
		for (i = 0; i < n; i++) {
			s2 += v(i,n);
		}
	}
	printf("  :  %.17g \t Diff: %g\n", s2, s2 - s1);
//	printf("  multi-threaded version            : %.17g \t Diff: %g \n", s2, s2 - s1);

	// SUMMATION USING REPRODUCIBLE TYPE
	dISetZero(I_s);
	for (i = 0; i < n; i++) {
		dIAddd(&I_s, v(i,n));
	}
	s1 = Iconv2d(I_s);
	printf("\n>>%20s :  %.17g \n", "Idouble [sequential]", s1);
	// multi-threaded
	dISetZero(I_s);
	#pragma omp parallel
	{
		I_double tmp;
		dISetZero(tmp);
		if (omp_get_thread_num() == 0)
			printf("            %2d theads", omp_get_num_threads());
		#pragma omp for schedule(dynamic, chunk) nowait
		for (i = 0; i < n; i++) {
			dIAddd(&tmp, v(i,n));
		}

		#pragma omp critical
		{
			dIAdd(&I_s, tmp);
		}
	}
	s2 = Iconv2d(I_s);
	printf("  :  %.17g \t Diff: %g\n", s2, s2 - s1);
}
