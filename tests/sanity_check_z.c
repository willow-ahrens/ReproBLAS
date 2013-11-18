#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "debug.h"

extern int sanity_check(int n, double* x, double* y, int* status);
extern const char* name();
extern const char* desc();

int main (int argc, char** args) {
	int status[6];
	int size = read_int(argc, args, "-n", 4096);
	int i;
	double* x = (double*)malloc(4 * size * sizeof(double));
	double* y = (double*)malloc(4 * size * sizeof(double));
	
	for (i = 0; i < 4*size; i++) {
		y[i] = 1.0;
	}

	fprintf(stdout, "|%8s | ", name());
	sanity_check(size, x, y, status);

	fprintf(stdout, "%s", desc());

	int left = 48 - strlen(desc());
	if (left > 0)
	for (i = 0; i < left; i++) fprintf(stdout, " ");

	int S = status[0] + status[1] + status[2] + status[3] + status[4];
	if (S == 0) {
		fprintf(stdout, " | PASSED |");
	}
	else {
		fprintf(stdout, " | FAILED |");
		if (status[0] != 0) {
			fprintf(stdout, "%8s", "SUM");
		}
		if (status[1] != 0) {
			fprintf(stdout, "%8s", "ASUM");
		}
		if (status[2] != 0) {
			fprintf(stdout, "%8s", "DOTC");
		}
		if (status[3] != 0) {
			fprintf(stdout, "%8s", "DOTU");
		}
		if (status[4] != 0) {
			fprintf(stdout, "%8s", "RDNRM2");
		}
	}
	fprintf(stdout, "\n");

	free(x);
	free(y);

	return 0;
}
