#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <IndexedFP.h>
#include "doubledouble.h"

double complex v(int i, int n) {
	return  sin(2*M_PI * (i / (double)n - 0.5))
		+ cos(2*M_PI * (i / (double)n - 0.5)) * _Complex_I ;
//	return 0.5 * i /n + 0.5 * _Complex_I;
}

int main (int argc, char** args) {
	int n = 1024;
	int i;
	double complex s1, s2;
	double ss[4];
	dIcomplex I_s;

	printf("\nSum of sin(PI * (i / n - 0.5) \n");
	// SUMMATION USING PRIMITIVE TYPE
	s1 = 0;
	for (i = 0; i < n; i++) {
		s1 += v(i,n);
	}
	printf("\n> Result with normal summation      : %.17g, %.17g \n", creal(s1), cimag(s1));
	// reverse order
	s2 = 0;
	for (i = n-1; i >= 0; i--) {
		s2 += v(i,n);
	}
	printf("  reversing order                   : %.17g, %.17g\n", creal(s2), cimag(s2));
	printf("  >> Diff                           : %.17g, %.17g\n", creal(s1 - s2), cimag(s1 - s2));

	// SUMMATION USING REPRODUCIBLE TYPE
	zISetZero(I_s);
	for (i = 0; i < n; i++) {
		zIAddz(&I_s, v(i,n));
	}
	s1 = Iconv2z(I_s);
	printf("\n> Result with reproducible summation: %.17g, %.17g \n", creal(s1), cimag(s1));
	// reverse order
	zISetZero(I_s);
	for (i = n-1; i >= 0; i--) {
		zIAddz(&I_s, v(i,n));
	}
	s2 = Iconv2z(I_s);
	printf("  reversing order                   : %.17g, %.17g \n", creal(s2), cimag(s2));
	printf("  >> Diff                           : %.17g, %.17g\n", creal(s1 - s2), cimag(s1 - s2));


}
