#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "debug.h"

extern int sanity_check_s(int n, float* x, float* y, int* status);
extern const char* name();
extern const char* desc();

int main (int argc, char** args) {
	int status[4];
	int size = read_int(argc, args, "-n", 4096);
	int i;
	float* x = (float*)malloc(size * sizeof(float));
	float* y = (float*)malloc(size * sizeof(float));
	
	for (i = 0; i < size; i++) {
		y[i] = 1.0;
	}

	fprintf(stdout, "|%8s | ", name());

	sanity_check_s(size, x, y, status);
	fprintf(stdout, "%s", desc());
	int left = 48 - strlen(desc());
	if (left > 0)
	for (i = 0; i < left; i++) fprintf(stdout, " ");
	int S = status[0] + status[1] + status[2] + status[3];
	if (S == 0) {
		fprintf(stdout, " | PASSED |");
	}
	else {
		fprintf(stdout, " | FAILED |");
//----------------------------------Possible misnaming? should be RSSUM, etc, etc?
		if (status[0] != 0) {
			fprintf(stdout, "%8s", "RDSUM");
		}
		if (status[1] != 0) {
			fprintf(stdout, "%8s", "RDASUM");
		}
		if (status[2] != 0) {
			fprintf(stdout, "%8s", "RDDOT");
		}
		if (status[3] != 0) {
			fprintf(stdout, "%8s", "RDNRM2");
		}
	}
	fprintf(stdout, "\n");

	free(x);
	free(y);

	return 0;
}
