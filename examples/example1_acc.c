#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <rblas.h>

static struct timeval start;
static struct timeval end;

void tic() {
  gettimeofday( &start, NULL );
}

double toc( )
{
  gettimeofday( &end, NULL );

  return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

#include <IndexedFP.h>

//double v(int i, int n) {
//	return  sin(2*M_PI * (i / (double)n - 0.5));
//}

#define v(i,n) sin(2 * M_PI * (i / (double)n - 0.5))

int main (int argc, char** args) {
	int n = 1000000;
	int i;
	double s1, s2;
	double ss[4];
	double elapsed;
	double t;
	dIAccum I_s;

	printf("\nSum of sin(2* M_PI * (i / n - 0.5).  n = %d \n", n);
	// SUMMATION USING PRIMITIVE TYPE
	s1 = 0;
	tic();
	for (i = 0; i < n; i++) {
		s1 += v(i,n);
	}
	elapsed = toc();
	printf("\n>>%16s :  %.17g \t [%5.1es] \n", "native double", s1, elapsed);
	// reverse order
	s2 = 0;
	for (i = n-1; i >= 0; i--) {
		s2 += v(i,n);
	}
	printf("  %16s :  %.17g \t Diff: %g \n", "reverse order", s2, s2 - s1);

	// SUMMATION USING REPRODUCIBLE TYPE
	tic();
	dInitAccum(&I_s, 256);
	for (i = 0; i < n; i++) {
		t = v(i,n);
		dAccumulate(&I_s, t);
	}
	s1 = dExtractAcc(&I_s);
	elapsed = toc();
	printf("\n>>%16s :  %.17g \t [%5.1es] \n", "dIAccum", s1, elapsed);
	// reverse order
	dResetAcc(&I_s);
	for (i = n-1; i >= 0; i--) {
		t = v(i,n);
		dAccumulate(&I_s, t);
	}
	s2 = dExtractAcc(&I_s);
	dDestroyAcc(&I_s);
	printf("  %16s :  %.17g \t Diff: %g \n", "reverse order", s2, s2 - s1);


	// SUMMATION USING double-double
	int fold = 2;
	tic();
	for (i = 0; i < 4; i++) ss[i] = 0.0;
	for (i = 0; i < n; i++) {
		dndpd(fold, ss, v(i,n));
	}
	s1 = ss[0] + ss[1];
	elapsed = toc();
	printf("\n>>%16s :  %.17g \t [%5.1es] \n", "double-double", s1, elapsed);
	//reverse
	for (i = 0; i < 4; i++) ss[i] = 0.0;
	for (i = n-1; i >= 0; i--) {
		dndpd(fold, ss, v(i,n));
	}
	s2 = ss[0] + ss[1];
	printf("  %16s :  %.17g \t Diff: %g \n", "reverse order", s2, s2 - s1);
}
