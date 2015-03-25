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
#define OUT_RCSUM  0
#define OUT_RSCASUM 1
#define OUT_RCDOTC  2
#define OUT_RSCNRM2 3

#define OUT_CSUM  4
#define OUT_SCASUM 5
#define OUT_CDOTC  6
#define OUT_SCNRM2 7

extern void dndpd(int, double*, double);

int main( int argc, char **argv ) {
	float complex* v;
	float complex* y;
	int incv;
	int incy = 1;
	int i;
	int NB;
	int n;
	int step = 1;
	int start;
	int stop;
	float ref;
	float complex zref;
	int dtype = 1;
	float cond = 1e3;
	int shift = 0;
	int from_file;
	char fname[128];
	float* prev_result = NULL;
	float out_result [NB_ARGOUT];
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
	cond  = read_float(argc, argv, "-K", 1e4);
	print_header = find_option(argc, argv, "--print-header") >= 0;

	if (find_option(argc, argv, "-K") >= 0)
		dtype = 5;

	start = stop = 128;
	step = 1;
	read_irange(argc, argv, "-n", &start, &stop, &step);

	// ALLOCATE MEMORY
	ref = 0.0;
	v = (float complex*)malloc(stop * sizeof(float complex));
	y = (float complex*)malloc(stop * sizeof(float complex));
//	sgenvec(stop, y, 0, 1.0);
	float* yptr = (float*) y;
	for (i = 0; i < stop; i++) yptr[i*2*incy] = yptr[i*2*incy+1] = 1.0;
	
	incv  = 1;
	NB    = read_int(argc, argv, "-nb", 1024);

	double ddR[6];
	double ddI[6];
	float zerr[2];
	float tmp_err[2];
	float factor = pow(2.0, shift);

	if (find_option(argc, argv, "-K") >= 0)
		dtype = 5;

	if (print_header) {
	printf("%8s", "N");
	#ifdef CALL_SCASUM
	printf("%10s", "SCASUM");
	#endif
	printf("%10s", "RSCASUM");
	printf("%10s", "RCSUM");
	#ifdef CALL_SCNRM2
	printf("%10s", "SCNRM2");
	#endif
	printf("%10s", "RSCNRM2");
	#ifdef CALL_CDOTU
	printf("%10s", "CDOTU");
	#endif
	printf("%10s", "RCDOTU");
	#ifdef CALL_CDOTC
	printf("%10s", "CDOTC");
	#endif
	printf("%10s", "RCDOTC");

	printf(" %10s", "RCOND(x)");
	printf(" %10s", "RCOND(x*y)");
	printf("  x[i] \n");
	}
	float ref_err, ref_rerr;
	float relerr;
	float rcond;
	int s;
	float* vptr = (float*) v;

	drandomseed();	

#define ACC_CHECK_	\
		ref_err = fabs((((ddR[0] - ref) + ddR[1]) + ddR[2]) + ddR[3]);	\
		ref_rerr = fabs(ddR[0])>fabs(ref)?fabs(ddR[0]):fabs(ref);	\
		if (ref_rerr != 0.0)	\
			ref_rerr = ref_err / ref_rerr;

#define QD_RELERR(DD,R,E) { \
		float abserr = fabs((((DD[0] - R) + DD[1]) + DD[2]) + DD[3]);	\
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
			sgenvec_(n, vptr, 2*incv, dtype, factor);
			sgenvec_(n, vptr + 1, 2*incv, dtype, factor);
		}

		printf("%8d", n);

		//==== SCASUM ====
		for (i = 0; i < 6; i++)
			ddR[i] = 0.0;

		for (i = 0; i < n; i++) {
			dndpd(4, ddR, fabs(vptr[i*incv*2]));
			dndpd(4, ddR, fabs(vptr[i*incv*2+1]));
		}
		rcond = ddR[0];
		
#		ifdef CALL_SCASUM
		CALL_SCASUM(ref, n, v, incv);
		out_result[OUT_SCASUM] = ref;
		ACC_CHECK
#		endif

		relerr = 0.0;
		for (s = 0; s < 53; s++) {
			ref = rscasum(n, v, incv);

			// log output
			if (s == 0)
				out_result[OUT_RSCASUM] = ref;

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


		//==== RCSUM ====
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
			zref = rcsum(n, v, incv);

			// log output
//			if (s == 0)
//				out_result[OUT_RCSUM] = ref;

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


		//==== RSCNRM2 ====
		for (i = 0; i < 6; i++)
			ddR[i] = 0.0;
		float scale = 0.0;
		float nscale = 0.0;
		int escale;
		int j;
		double ph, pl;
		double dnrm2[4];

		for (i = 0; i < n; i++) {
			float tmp = fabs(vptr[i*2*incv]);
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

		//==== SCNRM2 ====
#		ifdef CALL_SCNRM2
		CALL_SCNRM2(ref, n, v, incv);
		out_result[OUT_SCNRM2] = ref;

		ref = ref / scale;
		ph = TwoProd(ref, ref, &pl);
		for (i = 0; i < 4; i++) dnrm2[i] = ddR[i];
		dndpd(4, dnrm2, -ph);
		dndpd(4, dnrm2, -pl);
		ref_rerr = 2 * fabs(dnrm2[0]) / ddR[0];

		fprintf(stdout, "%10.2g", ref_rerr);	
#		endif

		ref = rscnrm2(n, v, incv);
		out_result[OUT_RSCNRM2] = ref;
		ref = ref / scale;
		ph = TwoProd(ref, ref, &pl);
		for (i = 0; i < 4; i++) dnrm2[i] = ddR[i];
		dndpd(4, dnrm2, -ph);
		dndpd(4, dnrm2, -pl);
		ref_rerr = 2 * fabs(dnrm2[0]) / ddR[0];

		fprintf(stdout, "%10.2g", ref_rerr);

		//==== CDOTU ====
		for (i = 0; i < 6; i++)
			ddR[i] = ddI[i] = 0.0;
		float rcondp = 0.0;
		float* yptr = (float*) y;
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


#		ifdef CALL_CDOTU
		CALL_CDOTU(zref, n, v, incv, y, incy);
		QD_RELERR(ddR,creal(zref), zerr[0]);
		QD_RELERR(ddI,cimag(zref), zerr[1]);
		fprintf(stdout, "%10.1e", zerr[0] < zerr[1] ? zerr[1] : zerr[0]);
#		endif

		zerr[0] = zerr[1] = 0.0;
		for (s = 0; s < 53; s++) {
			zref = rcdotu(n, v, incv, y, incy);

			// log output
			if (s == 0)
				out_result[OUT_RCDOTC] = ref;

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

		//==== CDOTC ====
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


#		ifdef CALL_CDOTC
		CALL_CDOTC(zref, n, v, incv, y, incy);
		QD_RELERR(ddR,creal(zref), zerr[0]);
		QD_RELERR(ddI,cimag(zref), zerr[1]);
		fprintf(stdout, "%10.1e", zerr[0] < zerr[1] ? zerr[1] : zerr[0]);
#		endif

		zerr[0] = zerr[1] = 0.0;
		for (s = 0; s < 53; s++) {
			zref = rcdotc(n, v, incv, y, incy);

			// log output
			if (s == 0)
				out_result[OUT_RCDOTC] = ref;

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
