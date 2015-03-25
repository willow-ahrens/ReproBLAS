#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <fenv.h>

#include <rblas.h>

#include "tmp_config.h"

#include "debug.h"

#define NB_ARGOUT  8 
#define OUT_RDSUM  0
#define OUT_RDASUM 1
#define OUT_RDDOT  2
#define OUT_RDNRM2 3

#define OUT_DSUM  4
#define OUT_DASUM 5
#define OUT_DDOT  6
#define OUT_DNRM2 7

extern void dndpd(int, double*, double);

int main( int argc, char **argv ) {
	double* v;
	double* y;
	int incv;
	int i;
	int NB;
	int n;
	int step = 1;
	int start;
	int stop;
	double ref;
	int IONE = 1;
	int dtype = 1;
	double cond = 1e3;
	int shift = 0;
	int from_file;
	char fname[128];
	double* prev_result = NULL;
	double out_result [NB_ARGOUT];
	int    check_prev_result = -1;
	int    write_output      = -1;
	int    print_header;

	if (find_option(argc, argv, "--help") >= 0) {
		printf("  OPTIONS:\n");
		show_opt("-n", "vector size");
		show_opt("-d", "input generating data type");
		show_opt("-K", "condition number of input vector");
		show_opt("-s", "shifting input data by a number of bits");
		show_opt("-f", "file name to read data");
		show_opt("-c", "Reference value to be checked with.");
		show_opt("-o", "Output file.");
		return 0;
	}

	from_file    = find_option(argc, argv, "-f");
	print_header = find_option(argc, argv, "--print-header") > 0;
	if (from_file >= 0 && from_file < argc - 1) {
		strcpy(fname, argv[from_file + 1]);
		from_file = 1;
		if (!fileExists(fname)) {
			fprintf(stderr, "FAILED: file %s not exists.\n", fname);
			return -1;
		}
		readFile(fname, &stop, (void**)&v, sizeof(double));
		start = stop;
		step = 1;
		y = (double*)malloc(stop * sizeof(double));
		for (i = 0; i < stop; i++) y[i] = 1.0;

		check_prev_result = find_option(argc, argv, "-c");
		if (check_prev_result >= argc - 1)
			check_prev_result = -1;

		if (check_prev_result >= 0) {
			strcpy(fname, argv[check_prev_result + 1]);
			// read

			if (!fileExists(fname)) {
				fprintf(stderr, "WARNING: file %s not exists.", fname);
				fprintf(stderr, " IGNORE reference input.\n");
				check_prev_result = -1;
			}
			else {
				int tmp;
				readFile(fname, &tmp, (void*)&prev_result, sizeof(double));
			}
		}


		write_output = find_option(argc, argv, "-o");
		if (write_output >= argc - 1)
			write_output = -1;

		if (write_output >= 0) {
			strcpy(fname, argv[write_output + 1]);
		}
	}
	else {
		from_file = 0;
		dtype = read_int(argc, argv, "-d", 1);
		shift = read_int(argc, argv, "-s", 0);
		cond  = read_double(argc, argv, "-K", 1e4);
		if (find_option(argc, argv, "-K") >= 0)
			dtype = 5;

		start = stop = 128;
		step = 1;
		read_irange(argc, argv, "-n", &start, &stop, &step);

		// ALLOCATE MEMORY
		ref = 0.0;
		v = (double*)malloc(stop * sizeof(double));
		y = (double*)malloc(stop * sizeof(double));
//		dgenvec(stop, y, 0, 1.0);
		for (i = 0; i < stop; i++) y[i] = 1.0;
	}
	
	incv  = 1;
	NB    = read_int(argc, argv, "-nb", 1024);

	double dd[6];
	double factor = pow(2.0, shift);

	if (find_option(argc, argv, "-K") >= 0)
		dtype = 5;

	if (print_header) {
	printf("%8s", "N");
	#ifdef CALL_DASUM
	printf("%10s", "DASUM");
	#endif
	printf("%10s", "RDASUM");
	printf("%10s", "RDSUM");
	#ifdef CALL_DNRM2
	printf("%10s", "DNRM2");
	#endif
	printf("%10s", "RDNRM2");
	#ifdef CALL_DDOT
	printf("%10s", "DDOT");
	#endif
	printf("%10s", "RDDOT");

	printf(" %10s", "RCOND(x)");
	printf(" %10s", "RCOND(x*y)");
	printf("  x[i] \n");
	}
	double ref_err, ref_rerr;
	double relerr;
	double rcond;
	int s;

	drandomseed();	

#define ACC_CHECK_	\
		ref_err = fabs((((dd[0] - ref) + dd[1]) + dd[2]) + dd[3]);	\
		ref_rerr = fabs(dd[0])>fabs(ref)?fabs(dd[0]):fabs(ref);	\
		if (ref_rerr != 0.0)	\
			ref_rerr = ref_err / ref_rerr;

#define ACC_CHECK	\
		ref_err = fabs((((dd[0] - ref) + dd[1]) + dd[2]) + dd[3]);	\
		ref_rerr = fabs(dd[0])>fabs(ref)?fabs(dd[0]):fabs(ref);	\
		if (ref_rerr != 0.0)	\
			ref_rerr = ref_err / ref_rerr;	\
		fprintf(stdout, "%10.2g", ref_rerr);	

	for (n = start; n < stop + step; n += step) {
		if (!from_file) {
			v[0] = cond;
			dgenvec(n, v, dtype, factor);
		}

		printf("%8d", n);

		//==== DASUM ====
		for (i = 0; i < 6; i++)
			dd[i] = 0.0;

		for (i = 0; i < n; i++) {
			dndpd(4, dd, fabs(v[i]));
		}
		rcond = dd[0];
		
#		ifdef CALL_DASUM
		CALL_DASUM(ref, n, v, IONE);
		out_result[OUT_DASUM] = ref;
		ACC_CHECK
#		endif

		relerr = 0.0;
		for (s = 0; s < 53; s++) {
			ref = rdasum(n, v, IONE);

			// log output
			if (s == 0)
				out_result[OUT_RDASUM] = ref;

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


		//==== RDSUM ====
		for (i = 0; i < 6; i++)
			dd[i] = 0.0;

		for (i = 0; i < n; i++) {
			dndpd(4, dd, v[i]);
		}
		if (rcond > 0) {
			rcond = fabs(dd[0]) / rcond;
		}

		relerr = 0.0;
		for (s = 0; s < 53; s++) {
			ref = rdsum(n, v, IONE);

			// log output
			if (s == 0)
				out_result[OUT_RDSUM] = ref;

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


		//==== RDNRM2 ====
		for (i = 0; i < 6; i++)
			dd[i] = 0.0;
		double scale = 0.0;
		double nscale = 0.0;
		int escale;
		int j;
		double ph, pl;
		double dnrm2[4];

		for (i = 0; i < n; i++) {
			if (fabs(v[i]) > 2*scale) {
				frexp(fabs(v[i]), &escale);
				nscale = ldexp(0.5, escale);
				scale  = scale / nscale;
				scale  = scale * scale;
				for (j = 0; j < 6; j++)
					dd[j] *= scale;
				scale = nscale;
			}
			double tmp = v[i];
			if (scale > 0)
				tmp = tmp / scale;
			ph = TwoProd(tmp, tmp, &pl);
			dndpd(4, dd, ph);
			dndpd(4, dd, pl);
		}

		//==== DNRM2 ====
#		ifdef CALL_DNRM2
		CALL_DNRM2(ref, n, v, IONE);
		out_result[OUT_DNRM2] = ref;

		ref = ref / scale;
		ph = TwoProd(ref, ref, &pl);
		for (i = 0; i < 4; i++) dnrm2[i] = dd[i];
		dndpd(4, dnrm2, -ph);
		dndpd(4, dnrm2, -pl);
		ref_rerr = 2 * fabs(dnrm2[0]) / dd[0];

		fprintf(stdout, "%10.2g", ref_rerr);	
#		endif

		ref = rdnrm2(n, v, IONE);
		out_result[OUT_RDNRM2] = ref;
		ref = ref / scale;
		ph = TwoProd(ref, ref, &pl);
		for (i = 0; i < 4; i++) dnrm2[i] = dd[i];
		dndpd(4, dnrm2, -ph);
		dndpd(4, dnrm2, -pl);
		ref_rerr = 2 * fabs(dnrm2[0]) / dd[0];

		fprintf(stdout, "%10.2g", ref_rerr);	
		//==== DDOT ====
		for (i = 0; i < 6; i++)
			dd[i] = 0.0;
		double rcondp = 0.0;

		for (i = 0; i < n; i++) {
			ph = TwoProd(v[i], y[i], &pl);
			rcondp += fabs(v[i] * y[i]);
			dndpd(4, dd, ph);
			dndpd(4, dd, pl);
		}
		if (rcondp != 0)
			rcondp = fabs(dd[0]) / rcondp;
//		printf("lo: %g\n", sl);


#		ifdef CALL_DDOT
		CALL_DDOT(ref, n, v, IONE, y, IONE);
		out_result[OUT_DDOT] = ref;
		ACC_CHECK;
#		endif

		relerr = 0.0;
		for (s = 0; s < 53; s++) {
			ref = rddot(n, v, IONE, y, IONE);

			// log output
			if (s == 0)
				out_result[OUT_RDDOT] = ref;

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

	if (check_prev_result >= 0) {
		printf("%8s", ".REP.");
		int count = 0;

#define CHECK_REP(IND) \
		if (out_result[IND] == prev_result[IND])	\
			printf("%10s", "MATCHED");	\
		else{		\
			printf("%10s", "FAILED");	\
			count ++;	\
		}

		#ifdef CALL_DASUM
		CHECK_REP(OUT_DASUM);
		#endif

		CHECK_REP(OUT_RDASUM);
		CHECK_REP(OUT_RDSUM);

		#ifdef CALL_DNRM2
		CHECK_REP(OUT_DNRM2);
		#endif

		CHECK_REP(OUT_RDNRM2);

		#ifdef CALL_DDOT
		CHECK_REP(OUT_DDOT);
		#endif

		CHECK_REP(OUT_RDDOT);

		printf("\n");

//----
		if (count > 0) {
		printf("%8s", ".PREV.");

#undef CHECK_REP
#define CHECK_REP(IND) \
		printf("  %8g", prev_result[IND]);

		#ifdef CALL_DASUM
		CHECK_REP(OUT_DASUM);
		#endif

		CHECK_REP(OUT_RDASUM);
		CHECK_REP(OUT_RDSUM);

		#ifdef CALL_DNRM2
		CHECK_REP(OUT_DNRM2);
		#endif

		CHECK_REP(OUT_RDNRM2);

		#ifdef CALL_DDOT
		CHECK_REP(OUT_DDOT);
		#endif

		CHECK_REP(OUT_RDDOT);

		printf("\n");
//----
		printf("%8s", ".DIFF.");

#undef CHECK_REP
#define CHECK_REP(IND) \
		printf("  %8.2g", prev_result[IND] - out_result[IND]);

		#ifdef CALL_DASUM
		CHECK_REP(OUT_DASUM);
		#endif

		CHECK_REP(OUT_RDASUM);
		CHECK_REP(OUT_RDSUM);

		#ifdef CALL_DNRM2
		CHECK_REP(OUT_DNRM2);
		#endif

		CHECK_REP(OUT_RDNRM2);

		#ifdef CALL_DDOT
		CHECK_REP(OUT_DDOT);
		#endif

		CHECK_REP(OUT_RDDOT);

		printf("\n");
		}
//----
		free(prev_result);
	}

	if (write_output >= 0) {
		writeFile(fname, NB_ARGOUT, out_result, sizeof(double));
	}

	free(v);
	free(y);
	return 0;
}
