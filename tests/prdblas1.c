#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <fenv.h>

#include <prblas.h>
#include <rblas.h>
#include "debug.h"

void MPI_NRM2(void *invec, void* inoutvec, int* len, MPI_Datatype *dtype) {
	int i;
	double* pin = (double*) invec;
	double* pout = (double*) inoutvec;
	for (i = 0; i < *len; i++) {
		pout[i] = sqrt(pout[i] * pout[i] + pin[i] * pin[i]);
	}
}

int main( int argc, char **argv ) {
	int rank, nprocs;
	int n;
	double* v;
	double* y;
	int incv;
	double sum;
	int i;
	int iters;
	double tic, toc;
	double FLOPs;
	double* result;
	double* elapsed;
	double* blas_result;
	double* blas_elapsed;
	int tind;
	int fold;
	int check, strong;
	double  refsum;
	int dtype;
	double K;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	n     = read_int(argc, argv, "-n", 1024);
	fold  = read_int(argc, argv, "-f", 3);
	strong = find_option(argc, argv, "--strong");
	check = find_option(argc, argv, "--check") >= 0;
	iters = read_int(argc, argv, "-i", 100);

	dtype = read_int(argc, argv, "-d", 1);
	K     = read_double(argc, argv, "-K", 1e4);
	if (find_option(argc, argv, "-K") >= 0)
		dtype = 5;
	incv  = 1;

	int N, first;
	if (strong < 0) {
		// WEAK SCALING
		N = n * nprocs;
		first = n * rank;
	}
	else {
		N = n;
		int q = N / nprocs;
		int r = N % nprocs;

		if (rank < r) {
			first = rank * (q + 1);
			n = q + 1;
		}
		else {
			first = r * (q + 1) + (rank - r) * q;
			n = q;
		}
	}

	// GENERATE DATA
	v     = (double*)malloc(n * sizeof(double));
	y     = (double*)malloc(n * sizeof(double));
	v[0]  = K;
	dgenvec(n, v, dtype, 1.0);
	for (i = 0; i < n; i++) y[i] = 1.0;

	elapsed = calloc(16 , sizeof(double));
	blas_elapsed = calloc(16 , sizeof(double));
	result = calloc(16 , sizeof(double));
	blas_result = calloc(16 , sizeof(double));

	if (rank == 0) {
		printf("%5s", "P");
		printf("%12s", "N");
		printf("%10s", "PRDDOT_R");
		printf("%10s", "PRDDOT_A");
		printf("%10s", "PRDASUM_R");
		printf("%10s", "PRDASUM_A");
		printf("%10s", "PRDSUM_R");
		printf("%10s", "PRDSUM_A");
		printf("%10s", "PRDNRM2_R");
		printf("%10s", "PRDNRM2_A");
		#ifdef CALL_DDOT
		printf("%10s", "DDOT_R");
		printf("%10s", "DDOT_A");
		#endif
		#ifdef CALL_DASUM
		printf("%10s", "DASUM_R");
		printf("%10s", "DASUM_A");
		#endif
		#ifdef CALL_DNRM2
		printf("%10s", "DNRM2_R");
		printf("%10s", "DNRM2_A");
		#endif
//		printf("P N DDOT_RED DDOT_ALL DASUM_RED DASUM_ALL DNRM2 DNR2M_ALL DSUM DSUM_ALL");
	}

	// EFFECTIVE FLOP COUNTS: N
	FLOPs = n;
	double MHz = iters * (FLOPs / 1e6);
	
	double* tsum;

	tsum = (double*) malloc((1 + 2 * fold) * sizeof(double));


#ifdef CALL_DDOT
	int ione = 1;
	double lsum;
	MPI_Op mpi_nrm2;
	MPI_Op mpi_sum;
	MPI_Op_create(MPI_NRM2, 1, &mpi_nrm2);
	if (rank == 0)
	printf(" BLAS_DDOT BLAS_DASUM BLAS_DNRM");
#endif
	if (rank == 0)
	printf("\n");

	for (i = 0; i < iters; i++) {
		// REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prddot(MPI_COMM_WORLD, 0, n, v, 1, y, 1);
		toc = read_timer();
			
		elapsed[0] += (toc - tic);
		result [0]  = sum;

		// ALL REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prddot(MPI_COMM_WORLD, -1, n, v, 1, y, 1);
		toc = read_timer();
			
		elapsed[1] += (toc - tic);
		result [1]  = sum;

#ifdef CALL_DDOT
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// local sum
		lsum = CALL_DDOT(n, v, ione, y, ione);

		//REDUCTION
		MPI_Reduce(&lsum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[0] += (toc - tic);
		blas_result [0]  = sum;
#endif

//=====
		// REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prdasum(MPI_COMM_WORLD, 0, n, v, 1);
		toc = read_timer();
			
		elapsed[2] += (toc - tic);
		result [2]  = sum;

		// ALL REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prdasum(MPI_COMM_WORLD, -1, n, v, 1);
		toc = read_timer();
			
		elapsed[3] += (toc - tic);
		result [3]  = sum;

#ifdef CALL_DDOT
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// local sum
		lsum = CALL_DASUM(n, v, ione);

		//REDUCTION
		MPI_Reduce(&lsum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[1] += (toc - tic);
		blas_result [1]  = sum;
#endif

//=====

		// REDUCE
		MPI_Barrier(MPI_COMM_WORLD);
		tic = read_timer();
		sum = prdnrm2(MPI_COMM_WORLD, 0, n, v, 1);
		toc = read_timer();
			
		elapsed[4] += (toc - tic);
		result [4]  = sum;

		// ALL REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prdnrm2(MPI_COMM_WORLD, -1, n, v, 1);
		toc = read_timer();
			
		elapsed[5] += (toc - tic);
		result [5]  = sum;

#ifdef CALL_DDOT
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// local sum
		lsum = CALL_DNRM2(n, v, ione);

		//REDUCTION
		MPI_Reduce(&lsum, &sum, 1, MPI_DOUBLE, mpi_nrm2, 0, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[2] += (toc - tic);
		blas_result [2]  = sum;
#endif

//====

		// REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prdsum(MPI_COMM_WORLD, 0, n, v, 1);
		toc = read_timer();
			
		elapsed[6] += (toc - tic);
		result [6]  = sum;

		// ALL REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prdsum(MPI_COMM_WORLD, -1, n, v, 1);
		toc = read_timer();
			
		elapsed[7] += (toc - tic);
		result [7]  = sum;

	}

	MPI_Barrier(MPI_COMM_WORLD);

	// PRINT OUT
	if (rank == 0) {
		printf("%5d", nprocs);
		printf("%12d", N);
		for (i = 0; i < 8; i++) {
			printf("%10.2g" , MHz / elapsed[i]);
		}
#ifdef CALL_DDOT
		for (i = 0; i < 3; i++) {
			printf("%10.2g" , MHz / blas_elapsed[i]);
		}
#endif
		printf("\n");
	}

	if (check) {
		double* lv = (double*) malloc((int)N * sizeof(double));
		double* ly = (double*) malloc((int)N * sizeof(double));
		double err, rerr;
		if (rank == 0)
			printf("%17s", "Reproducibility");

		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Allgather(v, n, MPI_DOUBLE, lv, n, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgather(y, n, MPI_DOUBLE, ly, n, MPI_DOUBLE, MPI_COMM_WORLD);

//#ifdef CALL_PRDDOT
		sum = rddot(N, lv, 1, ly, 1);
		err = RELERR(sum, result[1]);
		MPI_Reduce(&err, &rerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			printf("%10.2g", RELERR(sum, result[0]));
			printf("%10.2g", rerr);
		}		
//#endif

//#ifdef CALL_PRDASUM
		sum = rdasum(N, lv, 1);

		err = RELERR(sum, result[3]);
		MPI_Reduce(&err, &rerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			printf("%10.2g", RELERR(sum, result[2]));
			printf("%10.2g", rerr);
		}		
//#endif

//#ifdef CALL_PRDSUM
		sum = rdsum(N, lv, 1);

		err = RELERR(sum, result[7]);
		MPI_Reduce(&err, &rerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			printf("%10.2g", RELERR(sum, result[6]));
			printf("%10.2g", rerr);
		}		
//#endif

//#ifdef CALL_PRDNRM2
		sum = rdnrm2(N, lv, 1);
		err = RELERR(sum, result[5]);
		MPI_Reduce(&err, &rerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			printf("%10.2g", RELERR(sum, result[4]));
			printf("%10.2g", rerr);
		}		
//#endif

//=====
#ifdef CALL_DDOT

		sum = CALL_DDOT(N, lv, ione, ly, ione);
		err = RELERR(sum, blas_result[1]);
		MPI_Reduce(&err, &rerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		if (rank == 0) {
			printf("%10.2g", RELERR(sum, blas_result[0]));
			printf("%10.2g", rerr);
		}
#endif

#ifdef CALL_DASUM
		sum = CALL_DASUM(N, lv, ione);
		err = RELERR(sum, blas_result[3]);
		MPI_Reduce(&err, &rerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		if (rank == 0) {
			printf("%10.2g", RELERR(sum, blas_result[2]));
			printf("%10.2g", rerr);
		}
#endif

#ifdef CALL_DNRM2
		sum = CALL_DNRM2(N, lv, ione);

		err = RELERR(sum, blas_result[5]);
		MPI_Reduce(&err, &rerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		if (rank == 0) {
			printf("%10.2g", RELERR(sum, blas_result[4]));
			printf("%10.2g", rerr);
		}
#endif
		free(lv);
		free(ly);
		if (rank == 0)
		printf("\n");
	}

#ifdef CALL_DDOT
	MPI_Op_free(&mpi_nrm2);
#endif
	// FREE MEMORY
	free(v);
	free(blas_result);
	free(result);
	free(elapsed);
	free(blas_elapsed);
	free(tsum);
	MPI_Finalize();
	return 0;
}
