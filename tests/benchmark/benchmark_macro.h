#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <fenv.h>

#include "tmp_config.h"

#include "debug.h"

#define time_(R,MASK,BIT,FCT,N,V,INCV,FLOPs)	\
	if (MASK & (1 << BIT)) {	\
		int iters = 0;			\
		double tic, toc; 		\
		double elapsed	= 0.0;	\
\
		while (elapsed < 0.1 || iters < 1000) {	\
			iters++;	\
			tic = read_timer();	\
			R = FCT(N, V, INCV);	\
			toc   = read_timer();	\
			elapsed += toc - tic;	\
		}	\
		fprintf(stdout, "%10.1f",iters * (((FLOPs) / 1e6) / elapsed));	\
	}

#define time1_(R,MASK,BIT,FCT,N,V,INCV,FLOPs)	\
	if (MASK & (1 << BIT)) {	\
		int iters = 0;			\
		double tic, toc; 		\
		double elapsed	= 0.0;	\
\
		while (elapsed < 0.1 || iters < 1000) {	\
			iters++;	\
			tic = read_timer();	\
			FCT(R, N, V, INCV);	\
			toc   = read_timer();	\
			elapsed += toc - tic;	\
		}	\
		fprintf(stdout, "%10.1f",iters * (((FLOPs) / 1e6) / elapsed));	\
	}

#define time2_(R,MASK,BIT,FCT,N,X,INCX,Y,INCY,FLOPs)	\
	if (MASK & (1 << BIT)) {	\
		int iters = 0;			\
		double tic, toc; 		\
		int IONE = 1;			\
		double elapsed	= 0.0;	\
\
		while (elapsed < 0.1 || iters < 1000) {	\
			iters++;	\
			tic = read_timer();	\
			R = FCT(N, X, INCX, Y, INCY);	\
			toc   = read_timer();	\
			elapsed += toc - tic;	\
		}	\
		fprintf(stdout, "%10.1f",iters * (((FLOPs) / 1e6) / elapsed));	\
	}

#define time3_(R,MASK,BIT,FCT,N,X,INCX,Y,INCY,FLOPs)	\
	if (MASK & (1 << BIT)) {	\
		int iters = 0;			\
		double tic, toc; 		\
		int IONE = 1;			\
		double elapsed	= 0.0;	\
\
		while (elapsed < 0.1 || iters < 1000) {	\
			iters++;	\
			tic = read_timer();	\
			FCT(R, N, X, INCX, Y, INCY);	\
			toc   = read_timer();	\
			elapsed += toc - tic;	\
		}	\
		fprintf(stdout, "%10.1f",iters * (((FLOPs) / 1e6) / elapsed));	\
	}

#define CHECK_HEADER_() \
	int incv = 1, incy = 1;	\
	int step = 1, start, stop;	\
	int check;	\
	int dtype;	\
	int tests = 1023;	\
	int cmp_blas = 0;	\
	int print_header = 0;	\
\
	if (find_option(argc, argv, "--help") >= 0) {	\
		printf("  OPTIONS:\n");	\
		show_opt("-n", "vector size");	\
		show_opt("--incv", "vector stride (v)");	\
		show_opt("--incy", "vector stride (y)");	\
		show_opt("-d", "input generating data type");	\
		show_opt("--check", "Check accuracy and reproducibility");	\
		show_opt("--cmp", "Compare with BLAS if availble");	\
		show_opt("--header", "Printing header");	\
		return 0;	\
	}	\
\
	incv = read_int(argc, argv, "--incv", incv);	\
	incy = read_int(argc, argv, "--incy", incy);	\
	dtype = read_int(argc, argv, "-d", 1);	\
	check = find_option(argc, argv, "--check");	\
	cmp_blas = find_option(argc, argv, "--cmp") >= 0;	\
	print_header = find_option(argc, argv, "--header") >= 0;	\
\
	if (find_option(argc, argv, "--dot") >= 0) {	\
		if (tests == 1023) tests = 0;	\
		tests = tests | (1 << DOT_BIT);	\
	}	\
	if (find_option(argc, argv, "--rdot") >= 0) {	\
		if (tests == 1023) tests = 0;	\
		tests = tests | (1 << RDOT_BIT);	\
	}	\
	if (find_option(argc, argv, "--rsum") >= 0) {	\
		if (tests == 1023) tests = 0;	\
		tests = tests | (1 << RSUM_BIT);	\
	}	\
	if (find_option(argc, argv, "--asum") >= 0) {	\
		if (tests == 1023) tests = 0;	\
		tests = tests | (1 << ASUM_BIT);	\
	}	\
	if (find_option(argc, argv, "--rasum") >= 0) {	\
		if (tests == 1023) tests = 0;	\
		tests = tests | (1 << RASUM_BIT);	\
	}	\
	if (find_option(argc, argv, "--nrm2") >= 0) {	\
		if (tests == 1023) tests = 0;	\
		tests = tests | (1 << NRM2_BIT);	\
	}	\
	if (find_option(argc, argv, "--rnrm2") >= 0) {	\
		if (tests == 1023) tests = 0;	\
		tests = tests | (1 << RNRM2_BIT);	\
	}	\
\
	if (cmp_blas < 1) {	\
		tests = tests & ~(1 << ASUM_BIT);	\
		tests = tests & ~(1 << NRM2_BIT);	\
		tests = tests & ~(1 << DOT_BIT);	\
	}	\
\
	start = stop = 128;	\
	step = 1;	\
	read_irange(argc, argv, "-n", &start, &stop, &step);

