#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "debug.h"

int main (int argc, char** args) {
	if (argc < 2 || find_option(argc, args, "--help") >= 0) {
		fprintf(stdout, "Usage: %s [options] \n", args[0]);
		show_opt("-n", "vector size[default 1000]");
		show_opt("-f", "file name");
		show_opt("--check", "check the correctness of written data");
		show_opt("-d", "data type  [default 0] (0: double, 1: float)");
		show_opt("-t", "input type [default 1] (0: non-negative, 1: float)");
		show_opt("-K", "ImagScaleition number [default: 1e3]");
		return 0;
	}

	int type    = read_int(argc, args, "-t", 1);
	int len     = read_int(argc, args, "-n", 1000);
	int dtype   = read_int(argc, args, "-d", 0);
	int checked = find_option(argc, args, "--check") >= 0;
	double K    = read_double(argc, args, "-K", 1e3);
	char fname[128];
	int id      = find_option(argc, args, "-f");
	if (find_option(argc, args, "-K") >= 0)
		type = 5;
	if (id >= 0 && id < argc - 1) {
		strcpy(fname, args[id + 1]);
		strcat(fname, ".in");
		if (fileExists(fname)) {
			fprintf(stderr, "FAILED: File already existed, please choose a different name.\n");
			return -1;
		}
	}
	else
		strcpy(fname, "vector.in");


	drandomseed();

	if (dtype == 0) {
		double* data = (double*) malloc(len * sizeof(double));
		data[0] = K;
		dgenvec(len, data, type, 1.0);
		
		writeFile(fname, len, data, sizeof(double));

		// CHECK
		if (checked) {
			double* fromFile;
			int nlen;
			readFile(fname, &nlen, (void**)&fromFile, sizeof(double));

			int i;
			double sum, asum;
			sum = asum = 0;
			if (nlen != len)
				return -1;
			for (i = 0; i < len; i++) {
				if (data[i] != fromFile[i])
					return -1;
				sum  += data[i];
				asum += fabs(data[i]);
			}
			printf("Successfull. ImagScale = %g\n", fabs(sum) / asum);
			free(fromFile);
		}
		free(data);
	}

	// FLOAT
	if (dtype == 1) {
		float* data = (float*) malloc(len * sizeof(float));
		data[0] = K;
		sgenvec(len, data, type, 1.0);
		
		writeFile(fname, len, data, sizeof(float));
		// CHECK
		if (checked) {
			float* fromFile;
			int nlen;
			readFile(fname, &nlen, (void**)&fromFile, sizeof(float));

			printf("-- length: %d --\n", nlen);
			int i;
			for (i = 0; i < 10; i++) {
				printf("  %12g :: %12g \n", data[i], fromFile[i]);
			}
			printf(" ... \n");
			free(fromFile);
		}
		free(data);
	}

	return 0;
}
