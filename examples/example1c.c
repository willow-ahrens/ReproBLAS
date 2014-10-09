#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <IndexedFP.h>

float complex v(int i, int n) {
	return  sin(M_PI * (i / (float)n - 0.5))
		+ cos(M_PI * (i / (float)n - 0.5)) * _Complex_I ;
//	return 0.5 * i /n + 0.5 * _Complex_I;
}

int main (int argc, char** args) {
	int n = 1024;
	int i;
	float complex s1, s2;
	float ss[4];
	sIcomplex I_s;

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
	cISetZero(I_s);
	for (i = 0; i < n; i++) {
		cIAddc(&I_s, v(i,n));
	}
	s1 = Iconv2c(I_s);
	printf("\n> Result with reproducible summation: %.17g, %.17g \n", creal(s1), cimag(s1));
	// reverse order
	cISetZero(I_s);
	for (i = n-1; i >= 0; i--) {
		cIAddc(&I_s, v(i,n));
	}
	s2 = Iconv2c(I_s);
	printf("  reversing order                   : %.17g, %.17g \n", creal(s2), cimag(s2));
	printf("  >> Diff                           : %.17g, %.17g\n", creal(s1 - s2), cimag(s1 - s2));


}
