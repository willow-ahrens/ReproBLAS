#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <fenv.h>

#include <rblas.h>
#include "tmp_config.h"

#include "debug.h"

extern void dndpd(int, double*, double);

int main( int argc, char **argv ) {
	float* v;
	float* y;
	int incv;
	float sum;
	int i;
	float lsum[32];
	float tic, toc;
	float FLOPs, GHz;
	float elapsed;
	int NB;
	int n;
	int step = 1;
	int start;
	int stop;
	float ref;
	int IONE = 1;
	int check;
	int dtype;
	float cond;
	int shift;
	int    print_header;

	if (find_option(argc, argv, "--help") >= 0) {
		printf("  OPTIONS:\n");
		show_opt("-n", "vector size");
		show_opt("-d", "input generating data type");
		show_opt("-K", "condition number of input vector");
		show_opt("-s", "shifting input data by a number of bits");
		show_opt("--check", "Check accuracy and reproducibility");
		return 0;
	}

	dtype = read_int(argc, argv, "-d", 1);
	NB    = read_int(argc, argv, "-nb", 1024);
	shift = read_int(argc, argv, "-s", 0);
	cond  = read_float(argc, argv, "-K", 1e4);
	print_header = find_option(argc, argv, "--print-header") > 0;
	incv  = 1;

	if (find_option(argc, argv, "-K") >= 0)
		dtype = 5;

	start = stop = 128;
	step = 1;
	read_irange(argc, argv, "-n", &start, &stop, &step);

	// ALLOCATE MEMORY
	ref = 0.0;
	v = (float*)malloc(stop * sizeof(float));
	y = (float*)malloc(stop * sizeof(float));

//	sgenvec(stop, y, 0, 1.0);
	for (i = 0; i < stop; i++)
		y[i] = 1.0;

	double dd[6];
	float factor = pow(2.0, shift);

	if (print_header) {
	printf("%8s", "N");
	#ifdef CALL_SASUM
	printf("%10s", "SASUM");
	#endif
	#ifdef CALL_RSASUM
	printf("%10s", "RSASUM");
	#endif
	#ifdef CALL_RSSUM
	printf("%10s", "RSSUM");
	#endif
	#ifdef CALL_SNRM2
	printf("%10s", "SNRM2");
	#endif
	#ifdef CALL_RSNRM2
	printf("%10s", "RSNRM2");
	#endif
	#ifdef CALL_SDOT
	printf("%10s", "SDOT");
	#endif
	#ifdef CALL_RSDOT
	printf("%10s", "RSDOT");
	#endif

	printf(" %10s", "RCOND(x)");
	printf(" %10s", "RCOND(x*y)");
	printf("  x[i] \n");
	}
	float ref_err, ref_rerr;
	float relerr;
	double rcond;
	int s;

	drandomseed();	

#define ACC_CHECK_	\
		ref_err = fabs((dd[0] - ref) + dd[1]);					\
		ref_rerr = fabs(dd[0])>fabs(ref)?fabs(dd[0]):fabs(ref);	\
		if (ref_rerr != 0.0)	\
			ref_rerr = ref_err / ref_rerr;

#define ACC_CHECK	\
		ref_err = fabs((dd[0] - ref) + dd[1]);					\
		ref_rerr = fabs(dd[0])>fabs(ref)?fabs(dd[0]):fabs(ref);	\
		if (ref_rerr != 0.0)	\
			ref_rerr = ref_err / ref_rerr;	\
		fprintf(stdout, "%10.2g", ref_rerr);	

	for (n = start; n < stop + step; n += step) {
		v[0] = cond;
		sgenvec(n, v, dtype, factor);

		printf("%8d", n);
		FLOPs = n;

		//==== SASUM ====
		for (i = 0; i < 6; i++)
			dd[i] = 0.0;

		for (i = 0; i < n; i++) {
			dndpd(2, dd, (double)fabs(v[i]));
		}
		rcond = dd[0];
		
#		ifdef CALL_SASUM
		CALL_SASUM(ref, n, v, IONE);
		ACC_CHECK
#		endif

#		ifdef CALL_RSASUM
		ref = reproBLAS_rsasum(n, v, IONE);
		ACC_CHECK;
#		endif


		//==== RDSUM ====
		for (i = 0; i < 6; i++)
			dd[i] = 0.0;

		for (i = 0; i < n; i++) {
			dndpd(2, dd, v[i]);
		}
		if (rcond > 0) {
			rcond = fabs(dd[0]) / rcond;
		}

#		ifdef CALL_RSSUM
		relerr = 0.0;
		for (s = 0; s < 53; s++) {
			ref = reproBLAS_rssum(n, v, IONE);
			ACC_CHECK_;
			if (relerr < ref_rerr) relerr = ref_rerr;

			for (i = 0; i < 6; i++) dd[i] *= 2;
			for (i = 0; i < n; i++)
				v[i] *= 2.0;
		}
		fprintf(stdout, "%10.1e", relerr);

		relerr = pow(2.0, -53);
		for (i = 0; i < n; i++) {
			v[i] *= relerr;
		}
#		endif

		//==== DNRM2 ====
		for (i = 0; i < 6; i++)
			dd[i] = 0.0;

		for (i = 0; i < n; i++) {
			dndpd(2, dd, (double)v[i] * (double)v[i]);
		}
		dd[0] = sqrt(dd[0]);
		dd[1] = 0.0;

#		ifdef CALL_SNRM2
		CALL_SNRM2(ref, n, v, IONE);
		ACC_CHECK
#		endif

#		ifdef CALL_RSNRM2
		relerr = 0.0;
		for (s = 0; s < 53; s++) {
			ref = reproBLAS_rsnrm2(n, v, IONE);
			ACC_CHECK_;
			if (relerr < ref_rerr) relerr = ref_rerr;

			for (i = 0; i < 6; i++) dd[i] *= 2;
			for (i = 0; i < n; i++)
				v[i] *= 2.0;
		}
		fprintf(stdout, "%10.1e", relerr);

		relerr = pow(2.0, -53);
		for (i = 0; i < n; i++) {
			v[i] *= relerr;
		}
#		endif
		//==== SDOT ====
		for (i = 0; i < 6; i++)
			dd[i] = 0.0;

		double rcondp = 0.0;
		for (i = 0; i < n; i++) {
			dndpd(2, dd, (double)v[i] * (double)y[i]);
			rcondp += fabs(v[i] * y[i]);
		}

		if (rcondp > 0) {
			rcondp = fabs(dd[0]) / rcondp;
		}

#		ifdef CALL_SDOT
		CALL_SDOT(ref, n, v, IONE, y, IONE);
		ACC_CHECK;
#		endif

#		ifdef CALL_RSDOT
		relerr = 0.0;
		for (s = 0; s < 53; s++) {
			ref = reproBLAS_rsdot(n, v, IONE, y, IONE);
			ACC_CHECK_;
			if (relerr < ref_rerr) relerr = ref_rerr;

			for (i = 0; i < 6; i++) dd[i] *= 2;
			for (i = 0; i < n; i++)
				v[i] *= 2.0;
		}
		fprintf(stdout, "%10.1e", relerr);

		relerr = pow(2.0, -53);
		for (i = 0; i < n; i++) {
			v[i] *= relerr;
		}
#		endif
		/*
		ref = reproBLAS_rsdot(n, v, incv, y, incv);
		ACC_CHECK;
		*/
//		fprintf(stdout, "   %10g %10g %10g", ref, dd[0], dd[1]);

		fprintf(stdout, " %10.1e %8.1e", rcond, rcondp);

		if (dtype == 0) {
			fprintf(stdout, "    rand48()");
		}
		if (dtype == 1) {
			fprintf(stdout, "    2*rand48()-1");
		}
		if (dtype == 2) {
			fprintf(stdout, "    rand48()+rand48()-1");
		}
		if (dtype == 3) {
			fprintf(stdout, "    Normal distribution");
		}
		if (dtype == 4) {
			fprintf(stdout, "    sin(2*PI*i/N)");
		}

		fprintf(stdout, "\n");


	}

	free(v);
	free(y);
	return 0;
}
