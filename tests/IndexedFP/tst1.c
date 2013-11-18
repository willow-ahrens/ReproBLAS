#include <stdio.h>
#include <stdlib.h>
#include "../../src/IndexedFP/IndexedFP_Fast.h"
#include "../../src/RepBlas1.h"

int main(int argc, char** args) {
	fprintf(stdout, "Test 1 \n");
	double x[5];
	double ind[36];
	double y[5];
	double z[5];
	int i, j;
	int fold = 3;
	int n = 5;

	// Generate data
	x[0] = 1.1;
	x[1] = 1024.1;
	x[2] = -0.25 / 1024.0;
	x[3] = -1024.1;
	x[4] = 1.1 / (1024.0 * 1024.0 * 1024.0);

	for (i = 0; i < 5; i++) z[i] = 1024.0 * 1024.0 * 1024.0 * 1024.0;

	dconv2Indv_Fast(5, 1.0, x, 1, ind, 1, 6, fold, 0);
	Indconv2dv_Fast(fold, 5, ind, 1, 6, y, 1);

	for (i = 0; i < 5; i++) {
		fprintf(stdout, ": From %g to : \n", x[i]);
		fprintf(stdout, ":  ");
		for (j = 0; j < 2 * fold; j++) {
			fprintf(stdout, " %g", ind[i + j * 2 * fold]);
		}
		fprintf(stdout, "\n");
		fprintf(stdout, ":  back to %g [%g] \n", y[i], y[i] - x[i]);
	}

	fprintf(stdout, "-----------------------\n");
	dIndUpdatev_Fast(fold, 0, n, ind, 1, 6, z, 1);
	Indconv2dv_Fast(fold, 5, ind, 1, 6, y, 1);

	for (i = 0; i < 5; i++) {
		fprintf(stdout, ": Update by %g \n", z[i]);
		fprintf(stdout, ":  ");
		for (j = 0; j < 2 * fold; j++) {
			fprintf(stdout, " %g", ind[i + j * 2 * fold]);
		}
		fprintf(stdout, "\n");
		fprintf(stdout, ":  back to %g [%g] \n", y[i], y[i] - x[i]);
	}

	fprintf(stdout, "-----------------------\n");
	fprintf(stdout, " daxpIy \n");
	daxpIy(n, fold, 1.0, x, 1, ind, 1, 6);

	Indconv2dv_Fast(fold, 5, ind, 1, 6, y, 1);

	for (i = 0; i < 5; i++) {
		fprintf(stdout, ":  ");
		for (j = 0; j < 2 * fold; j++) {
			fprintf(stdout, " %g", ind[i + j * 2 * fold]);
		}
		fprintf(stdout, "\n");
		fprintf(stdout, ":  back to %g [%g] \n", y[i], y[i] - 2 * x[i]);
	}
}
