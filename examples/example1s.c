#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <IndexedFP.h>

float v(int i, int n) {
	return  sin(M_PI * (i / (float)n - 0.5));
}

int main (int argc, char** args) {
	int n = 1024;
	int i;
	float s1, s2;
	double ss[4];
	I_float I_s;

	printf("\nSum of sin(PI * (i / n - 0.5) \n");
	// SUMMATION USING PRIMITIVE TYPE
	s1 = 0;
	for (i = 0; i < n; i++) {
		s1 += v(i,n);
	}
	printf("\n> Result with normal summation      : %.17g \n", s1);
	// reverse order
	s2 = 0;
	for (i = n-1; i >= 0; i--) {
		s2 += v(i,n);
	}
	printf("  reversing order                   : %.17g \t Diff: %g \n", s2, s2 - s1);

	// SUMMATION USING REPRODUCIBLE TYPE
	dISetZero(I_s);
	for (i = 0; i < n; i++) {
		sIAddf(&I_s, v(i,n));
	}
	s1 = Iconv2f(I_s);
	printf("\n> Result with reproducible summation: %.17g \n", s1);
	// reverse order
	dISetZero(I_s);
	for (i = n-1; i >= 0; i--) {
		sIAddf(&I_s, v(i,n));
	}
	s2 = Iconv2f(I_s);
	printf("  reversing order                   : %.17g \t Diff: %g\n", s2, s2 - s1);


	// SUMMATION USING float-float
	int fold = 2;
	for (i = 0; i < 4; i++) ss[i] = 0.0;
	for (i = 0; i < n; i++) {
		dndpd(fold, ss, v(i,n));
	}
	s1 = ss[0] + ss[1];
	printf("\n> Result with higher precision      : %.17g \n", s1);
	//reverse
	for (i = 0; i < 4; i++) ss[i] = 0.0;
	for (i = n-1; i >= 0; i--) {
		dndpd(fold, ss, v(i,n));
	}
	s2 = ss[0] + ss[1];
	printf("  reversing order                   : %.17g \t Diff: %g \n", s2, s2 - s1);
}
