#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <fenv.h>

#include <rblas.h>

#include "../src/types.h"

#include "tmp_config.h"

#include "debug.h"

#define NB_ARGOUT  8 
#define OUT_RZSUM  0
#define OUT_RDZASUM 1
#define OUT_RZDOTC  2
#define OUT_RDZNRM2 3

#define OUT_ZSUM  4
#define OUT_DZASUM 5
#define OUT_ZDOTC  6
#define OUT_DZNRM2 7

extern void dndpd(int, double*, double);

int main( int argc, char **argv ) {
	double complex* v;
	double complex* y;
	int incv;
	int incy = 1;
	int i;
	int NB;
	int n;
	int step = 1;
	int start;
	int stop;
	double ref;
	double complex zref;
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

	from_file = 0;
	dtype = read_int(argc, argv, "-d", 1);
	shift = read_int(argc, argv, "-s", 0);
	cond  = read_double(argc, argv, "-K", 1e4);
	print_header = find_option(argc, argv, "--print-header") >= 0;

	if (find_option(argc, argv, "-K") >= 0)
		dtype = 5;

	start = stop = 128;
	step = 1;
	read_irange(argc, argv, "-n", &start, &stop, &step);

	// ALLOCATE MEMORY
	ref = 0.0;
	v = (double complex*)malloc(stop * sizeof(double complex));
	y = (double complex*)malloc(stop * sizeof(double complex));
//	dgenvec(stop, y, 0, 1.0);
	double* yptr = (double*) y;
//	for (i = 0; i < stop; i++) {
//		yptr[i*2*incy] = yptr[i*2*incy+1] = 1.0;
//	}
	dgenvec_(stop, yptr, 2 * incy, 0, 1.0);
	dgenvec_(stop, yptr + 1, 2 * incy, 0, 1.0);
	
	incv  = 1;
	NB    = read_int(argc, argv, "-nb", 1024);

	double ddR[6];
	double ddI[6];
	double zerr[2];
	double tmp_err[2];
	double factor = pow(2.0, shift);

	if (find_option(argc, argv, "-K") >= 0)
		dtype = 5;

	if (print_header) {
	printf("%8s", "N");
	#ifdef CALL_DZASUM
	printf("%10s", "DZASUM");
	#endif
	printf("%10s", "RDZASUM");
	printf("%10s", "RZSUM");
	#ifdef CALL_DZNRM2
	printf("%10s", "DZNRM2");
	#endif
	printf("%10s", "RDZNRM2");
	#ifdef CALL_ZDOTU
	printf("%10s", "ZDOTU");
	#endif
	printf("%10s", "RZDOTU");
	#ifdef CALL_ZDOTC
	printf("%10s", "ZDOTC");
	#endif
	printf("%10s", "RZDOTC");

	printf(" %10s", "RCOND(x)");
	printf(" %10s", "RCOND(x*y)");
	printf("  x[i] \n");
	}
	double ref_err, ref_rerr;
	double relerr;
	double rcond;
	int s;
	double* vptr = (double*) v;

	drandomseed();	

#define ACC_CHECK_	\
		ref_err = fabs((((ddR[0] - ref) + ddR[1]) + ddR[2]) + ddR[3]);	\
		ref_rerr = fabs(ddR[0])>fabs(ref)?fabs(ddR[0]):fabs(ref);	\
		if (ref_rerr != 0.0)	\
			ref_rerr = ref_err / ref_rerr;

#define QD_RELERR(DD,R,E) { \
		double abserr = fabs((((DD[0] - R) + DD[1]) + DD[2]) + DD[3]);	\
		E = fabs(DD[0])>fabs(R)?fabs(DD[0]):fabs(R);	\
		if (E != 0.0)	\
			E = abserr / E;	\
	}

#define ACC_CHECK	\
		ref_err = fabs((((ddR[0] - ref) + ddR[1]) + ddR[2]) + ddR[3]);	\
		ref_rerr = fabs(ddR[0])>fabs(ref)?fabs(ddR[0]):fabs(ref);	\
		if (ref_rerr != 0.0)	\
			ref_rerr = ref_err / ref_rerr;	\
		fprintf(stdout, "%10.2g", ref_rerr);	

	for (n = start; n < stop + step; n += step) {
		if (!from_file) {
			v[0] = cond;
			dgenvec_(n, vptr, 2*incv, dtype, factor);
			dgenvec_(n, vptr + 1, 2*incv, dtype, factor);
		}

		printf("%8d", n);

		//==== DZASUM ====
		for (i = 0; i < 6; i++)
			ddR[i] = 0.0;

		for (i = 0; i < n; i++) {
			dndpd(4, ddR, fabs(vptr[i*incv*2]));
			dndpd(4, ddR, fabs(vptr[i*incv*2+1]));
		}
		rcond = ddR[0];
		
#		ifdef CALL_DZASUM
		CALL_DZASUM(ref, n, v, incv);
		out_result[OUT_DZASUM] = ref;
		ACC_CHECK
#		endif

		relerr = 0.0;
		for (s = 0; s < 53; s++) {
			ref = reproBLAS_rdzasum(n, v, incv);

			// log output
			if (s == 0)
				out_result[OUT_RDZASUM] = ref;

			ACC_CHECK_;
			if (relerr < ref_rerr) relerr = ref_rerr;

			for (i = 0; i < 6; i++) ddR[i] *= 2;
			for (i = 0; i < n; i++)
				v[i*incv] *= 2.0;
		}
		fprintf(stdout, "%10.1e", relerr);

		relerr = pow(2.0, -53);
		for (i = 0; i < n; i++) {
			v[i*incv] *= relerr;
		}


		//==== RZSUM ====
		for (i = 0; i < 6; i++)
			ddR[i] = ddI[i] = 0.0;

		for (i = 0; i < n; i++) {
			dndpd(4, ddR, vptr[i*2*incv]);
			dndpd(4, ddI, vptr[i*2*incv+1]);
		}
		if (rcond > 0) {
			rcond = fabs(ddR[0]) / rcond;
		}

		zerr[0] = zerr[1] = 0.0;

		for (s = 0; s < 53; s++) {
			zref = rzsum(n, v, incv);

			// log output
//			if (s == 0)
//				out_result[OUT_RZSUM] = ref;

//			ACC_CHECK_;
			QD_RELERR(ddR,creal(zref), tmp_err[0]);
			QD_RELERR(ddI,cimag(zref), tmp_err[1]);

			if (zerr[0] < tmp_err[0]) zerr[0] = tmp_err[0];
			if (zerr[1] < tmp_err[1]) zerr[1] = tmp_err[1];

			for (i = 0; i < 6; i++) {
				ddR[i] *= 2;
				ddI[i] *= 2;
			}
			for (i = 0; i < n; i++)
				v[i*incv] *= 2.0;
		}
		fprintf(stdout, "%10.1e", zerr[0] < zerr[1] ? zerr[1] : zerr[0]);

		relerr = pow(2.0, -53);
		for (i = 0; i < n; i++) {
			v[i*incv] *= relerr;
		}


		//==== RDZNRM2 ====
		for (i = 0; i < 6; i++)
			ddR[i] = 0.0;
		double scale = 0.0;
		double nscale = 0.0;
		int escale;
		int j;
		double ph, pl;
		double dnrm2[4];

		for (i = 0; i < n; i++) {
			double tmp = fabs(vptr[i*2*incv]);
			if (fabs(vptr[i*2*incv+1]) > tmp)
				tmp = fabs(vptr[i*2*incv+1]);
			if (tmp > 2*scale) {
				frexp(tmp, &escale);
				nscale = ldexp(0.5, escale);
				scale  = scale / nscale;
				scale  = scale * scale;
				for (j = 0; j < 6; j++)
					ddR[j] *= scale;
				scale = nscale;
			}
			tmp = vptr[i*2*incv];
			if (scale > 0)
				tmp = tmp / scale;
			ph = TwoProd(tmp, tmp, &pl);
			dndpd(4, ddR, ph);
			dndpd(4, ddR, pl);

			tmp = vptr[i*2*incv + 1];
			if (scale > 0)
				tmp = tmp / scale;
			ph = TwoProd(tmp, tmp, &pl);
			dndpd(4, ddR, ph);
			dndpd(4, ddR, pl);
		}

		//==== DZNRM2 ====
#		ifdef CALL_DZNRM2
		CALL_DZNRM2(ref, n, v, incv);
		out_result[OUT_DZNRM2] = ref;

		ref = ref / scale;
		ph = TwoProd(ref, ref, &pl);
		for (i = 0; i < 4; i++) dnrm2[i] = ddR[i];
		dndpd(4, dnrm2, -ph);
		dndpd(4, dnrm2, -pl);
		ref_rerr = 2 * fabs(dnrm2[0]) / ddR[0];

		fprintf(stdout, "%10.2g", ref_rerr);	
#		endif

		ref = reproBLAS_rdznrm2(n, v, incv);
		out_result[OUT_RDZNRM2] = ref;
		ref = ref / scale;
		ph = TwoProd(ref, ref, &pl);
		for (i = 0; i < 4; i++) dnrm2[i] = ddR[i];
		dndpd(4, dnrm2, -ph);
		dndpd(4, dnrm2, -pl);
		ref_rerr = 2 * fabs(dnrm2[0]) / ddR[0];

		fprintf(stdout, "%10.2g", ref_rerr);

		//==== ZDOTU ====
		for (i = 0; i < 6; i++)
			ddR[i] = ddI[i] = 0.0;
		double rcondp = 0.0;
		for (i = 0; i < n; i++) {
			ph = TwoProd(vptr[i*2*incv], yptr[i*2*incy], &pl);
			rcondp += fabs(vptr[i*2*incv] * yptr[i*2*incy]);
			dndpd(4, ddR, ph);
			dndpd(4, ddR, pl);

			ph = TwoProd(vptr[i*2*incv+1], -yptr[i*2*incy+1], &pl);
			rcondp += fabs(v[i*2*incv] * yptr[i*2*incy]);
			dndpd(4, ddR, ph);
			dndpd(4, ddR, pl);

			//---
			ph = TwoProd(vptr[i*2*incv], yptr[i*2*incy+1], &pl);
			dndpd(4, ddI, ph);
			dndpd(4, ddI, pl);

			ph = TwoProd(vptr[i*2*incv+1], yptr[i*2*incy], &pl);
			dndpd(4, ddI, ph);
			dndpd(4, ddI, pl);

		}
		if (rcondp != 0)
			rcondp = fabs(ddR[0]) / rcondp;


#		ifdef CALL_ZDOTU
		CALL_ZDOTU(zref, n, v, incv, y, incy);
		QD_RELERR(ddR,creal(zref), zerr[0]);
		QD_RELERR(ddI,cimag(zref), zerr[1]);
		fprintf(stdout, "%10.1e", zerr[0] < zerr[1] ? zerr[1] : zerr[0]);
#		endif

		zerr[0] = zerr[1] = 0.0;
		for (s = 0; s < 53; s++) {
			zref = rzdotu(n, v, incv, y, incy);

			// log output
			if (s == 0)
				out_result[OUT_RZDOTC] = ref;

			QD_RELERR(ddR,creal(zref), tmp_err[0]);
			QD_RELERR(ddI,cimag(zref), tmp_err[1]);
			if (relerr < ref_rerr) relerr = ref_rerr;

			if (zerr[0] < tmp_err[0]) zerr[0] = tmp_err[0];
			if (zerr[1] < tmp_err[1]) zerr[1] = tmp_err[1];

			for (i = 0; i < 6; i++) {
				ddR[i] *= 2;
				ddI[i] *= 2;
			}
			for (i = 0; i < n; i++)
				v[i*incv] *= 2.0;
		}
		fprintf(stdout, "%10.1e", zerr[0] < zerr[1] ? zerr[1] : zerr[0]);

		relerr = pow(2.0, -53);
		for (i = 0; i < n; i++) {
			v[i*incv] *= relerr;
		}

		//==== ZDOTC ====
		for (i = 0; i < 6; i++)
			ddR[i] = ddI[i] = 0.0;
		rcondp = 0.0;
		for (i = 0; i < n; i++) {
			ph = TwoProd(vptr[i*2*incv], yptr[i*2*incy], &pl);
			rcondp += fabs(vptr[i*2*incv] * yptr[i*2*incy]);
			dndpd(4, ddR, ph);
			dndpd(4, ddR, pl);

			ph = TwoProd(vptr[i*2*incv+1], yptr[i*2*incy+1], &pl);
			rcondp += fabs(v[i*2*incv] * yptr[i*2*incy]);
			dndpd(4, ddR, ph);
			dndpd(4, ddR, pl);

			//---
			ph = TwoProd(vptr[i*2*incv], yptr[i*2*incy+1], &pl);
			dndpd(4, ddI, ph);
			dndpd(4, ddI, pl);

			ph = TwoProd(vptr[i*2*incv+1], -yptr[i*2*incy], &pl);
			dndpd(4, ddI, ph);
			dndpd(4, ddI, pl);

		}
		if (rcondp != 0)
			rcondp = fabs(ddR[0]) / rcondp;


#		ifdef CALL_ZDOTC
		CALL_ZDOTC(zref, n, v, incv, y, incy);
		QD_RELERR(ddR,creal(zref), zerr[0]);
		QD_RELERR(ddI,cimag(zref), zerr[1]);
		fprintf(stdout, "%10.1e", zerr[0] < zerr[1] ? zerr[1] : zerr[0]);
#		endif

		zerr[0] = zerr[1] = 0.0;
		for (s = 0; s < 53; s++) {
			zref = rzdotc(n, v, incv, y, incy);

			// log output
			if (s == 0)
				out_result[OUT_RZDOTC] = ref;

			QD_RELERR(ddR,creal(zref), tmp_err[0]);
			QD_RELERR(ddI,cimag(zref), tmp_err[1]);
			if (relerr < ref_rerr) relerr = ref_rerr;

			if (zerr[0] < tmp_err[0]) zerr[0] = tmp_err[0];
			if (zerr[1] < tmp_err[1]) zerr[1] = tmp_err[1];

			for (i = 0; i < 6; i++) {
				ddR[i] *= 2;
				ddI[i] *= 2;
			}
			for (i = 0; i < n; i++)
				v[i*incv] *= 2.0;
		}
//		fprintf(stdout, "%10.1e|%10.1e", zerr[0] , zerr[1]);
		fprintf(stdout, "%10.1e", zerr[0] < zerr[1] ? zerr[1] : zerr[0]);

		relerr = pow(2.0, -53);
		for (i = 0; i < n; i++) {
			v[i*incv] *= relerr;
		}

		//=====
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
