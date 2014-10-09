#include <IndexedFP.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char** args) {
	Idouble idouble;
	Idouble I2;
	Ifloat ifloat;

	double dd[4];
	float  ff[4];

	dd[0] = 0;
	dd[1] = 1.0;
	dd[2] = drand48();
	dd[3] = -dd[2];

	double dtmp;
	double x = drand48();

	idouble = dconv2I(x);
	dISet(I2, idouble);
	dINeg(I2);
	dIAdd(&I2, idouble);
	printf("x + (-x) = %g\n", Iconv2d(I2));

	int i;

	// conversion
	for (i = 0; i < 4; i++) {
		idouble = dconv2I(dd[i]);
		dtmp = Iconv2d(idouble);
		printf("x = %g -> dconv2I(x) = ", dd[i]);
		dIprint(idouble);
		printf("\nIconv2d(dconv2I(x)) = %g. err = %g \n", dtmp, dd[i] - dtmp);
	}

	// summing
	dISetZero(idouble);
	dtmp = 0.0;
	for (i = 0; i < 4; i++) {
		dIAddd(&idouble, dd[i]);
		dtmp += dd[i];
	}
	printf("normal sum: %g\n", dtmp);
	printf("reproducible sum: %g\n", Iconv2d(idouble));

	// summing backward
	dISetZero(idouble);
	dtmp = 0.0;
	for (i = 3; i >= 0; i--) {
		dIAddd(&idouble, dd[i]);
		dtmp += dd[i];
	}
	printf("normal sum: %g\n", dtmp);
	printf("reproducible sum: %g\n", Iconv2d(idouble));

	// ARRAY OF INDEXED
	Idouble idd[4];
	for (i = 0; i < 4; i++) {
		idd[i] = dconv2I(dd[i]);
	}
	// SUM
	dISetZero(idouble);
	dtmp = 0.0;
	for (i = 3; i >= 0; i--) {
		dIAdd(&idouble, idd[i]);
	}
	printf("reproducible sum: %g\n", Iconv2d(idouble));

}
