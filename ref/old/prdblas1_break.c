#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <fenv.h>

#include <prblas.h>
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
	int i, j;
	int iters;
	double tic, toc, toc1;
	double FLOPs;
	double* result;
	double* elapsed;
	double* blas_result;
	double* blas_elapsed;
	int tind;
	int fold;
	int check;
	double  refsum;
	int dtype;

	MPI_Init(&argc, &argv);
	RMPI_Init();
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	n     = read_int(argc, argv, "-n", 1024);
	iters = read_int(argc, argv, "-i", 100);
	fold  = read_int(argc, argv, "-f", 3);
	check = find_option(argc, argv, "--strong");
	dtype = read_int(argc, argv, "-d", 1);
	incv  = 1;

	int N, first;
	if (check < 0) {
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
	dgenvec(n, v, dtype, 1.0);
	for (i = 0; i < n; i++) y[i] = 1.0;

	if (rank == 0)
		printf("P N DDOT_RED DDOT_ALL DASUM_RED DASUM_ALL DNRM2 DNR2M_ALL DSUM DSUM_ALL");

	// EFFECTIVE FLOP COUNTS: N
	FLOPs = n;
	double MHz = iters * (FLOPs / 1e6);
	
	I_double tsum, gsum;

//	tsum = (double*) malloc((2 + 2 * fold) * sizeof(double));
//	gsum = (double*) malloc((2 + 2 * fold) * sizeof(double));

	elapsed = calloc(16 , sizeof(double));
	double* comm = calloc(16 , sizeof(double));
	double* blas_comm = calloc(16 , sizeof(double));
	blas_elapsed = calloc(16 , sizeof(double));
	result = calloc(16 , sizeof(double));
	blas_result = calloc(16 , sizeof(double));

	double lsum;
	int ione = 1;
#ifdef CALL_DDOT
	MPI_Op mpi_nrm2;
	MPI_Op_create(MPI_NRM2, 1, &mpi_nrm2);
#endif

#ifdef CALL_DASUM
	lsum = CALL_DASUM(n, v, ione);
#endif
#ifdef CALL_DNRM2
	lsum = CALL_DNRM2(n, v, ione);
#endif
#ifdef CALL_DDOT
	lsum = CALL_DDOT(n, v, ione, y, ione);
	if (rank == 0)
	printf(" BLAS_DDOT BLAS_DASUM BLAS_DNRM\n");
#endif
	if (rank == 0)
	printf("\n");


	prddot(MPI_COMM_WORLD, 0, n, v, 1, y, 1);
	prdasum(MPI_COMM_WORLD, 0, n, v, 1);
	prdsum(MPI_COMM_WORLD, 0, n, v, 1);
	prdnrm2(MPI_COMM_WORLD, 0, n, v, 1);

	for (i = 0; i < iters; i++) {
		// REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// LOCAL COMPUTATION
		tsum = ddotI(n, v, 1, y, 1);
		toc1 = read_timer();

		// REDUCTION
		MPI_Reduce(&tsum, &gsum, 1, MPI_IDOUBLE, MPI_RSUM, 0, MPI_COMM_WORLD);
		result[0]  = Iconv2d(gsum);
		toc = read_timer();
			
		elapsed[0] += (toc - tic);
		comm[0] += toc1 - tic;

#ifdef CALL_DDOT
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// local sum
		lsum = CALL_DDOT(n, v, ione, y, ione);

		//REDUCTION
		toc1 = read_timer();
		MPI_Reduce(&lsum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[0] += (toc - tic);
		blas_comm[0] += (toc - toc1);
		blas_result[0]  = sum;
#endif

//=====
		// REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// LOCAL COMPUTATION
		tsum = dasumI(n, v, 1);
		toc1 = read_timer();

		// REDUCTION
		MPI_Reduce(&tsum, &gsum, 1, MPI_IDOUBLE, MPI_RSUM, 0, MPI_COMM_WORLD);
		result[1]  = Iconv2d(gsum);
		toc = read_timer();
			
		elapsed[1] += (toc - tic);
		comm[1] += toc1 - tic;

#ifdef CALL_DDOT
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// local sum
		lsum = CALL_DASUM(n, v, ione);

		//REDUCTION
		toc1 = read_timer();
		MPI_Reduce(&lsum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[1] += (toc - tic);
		blas_comm[1] += (toc - toc1);
		blas_result[1]  = sum;
#endif

//=====

		// REDUCE
		
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// LOCAL COMPUTATION
		dISetZero(tsum);
		double scale = dnrm2I(n, v, 1, &tsum);
		toc1 = read_timer();

		// REDUCTION
		MPI_Reduce(&tsum, &gsum, 1, MPI_IDOUBLE_SCALE, MPI_RNRM2, 0, MPI_COMM_WORLD);
		result[2]  = Iconv2d(gsum);
		toc = read_timer();
			
		elapsed[2] += (toc - tic);
		comm[2] += toc1 - tic;

#ifdef CALL_DDOT
		MPI_Barrier(MPI_COMM_WORLD);
		lsum = CALL_DNRM2(n, v, ione);

		tic = read_timer();
		// local sum
		lsum = CALL_DNRM2(n, v, ione);

		//REDUCTION
		toc1 = read_timer();
		MPI_Reduce(&lsum, &sum, 1, MPI_DOUBLE, mpi_nrm2, 0, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[2] += (toc - tic);
		blas_comm[2] += (toc - toc1);
		blas_result[2]  = sum;
#endif

//====

		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// LOCAL COMPUTATION
		tsum = dsumI(n, v, 1);
		toc1 = read_timer();

		// REDUCTION
		MPI_Reduce(&tsum, &gsum, 1, MPI_IDOUBLE, MPI_RSUM, 0, MPI_COMM_WORLD);
		result[3]  = Iconv2d(gsum);
		toc = read_timer();
			
		elapsed[3] += (toc - tic);
		comm[3] += toc1 - tic;

	}

	MPI_Barrier(MPI_COMM_WORLD);

	// PRINT OUT
	if (rank == 0) {
		printf("%5d", nprocs);
		printf("%12d", N);
		for (i = 0; i < 4; i++) {
			printf("%8.2f %12.3g %12.3g " , MHz / elapsed[i], elapsed[i], comm[i]);
		}
#ifdef CALL_DDOT
		for (i = 0; i < 3; i++) {
			printf("%8.2f %12.3g %12.3g" , MHz / blas_elapsed[i], blas_elapsed[i], blas_comm[i]);
		}
#endif
		printf("\n");
	}
	
#ifdef CALL_DDOT
	MPI_Op_free(&mpi_nrm2);
#endif
	// FREE MEMORY
	free(v);
	free(blas_result);
	free(result);
	free(comm);
	free(blas_comm);
	free(elapsed);
	free(blas_elapsed);
	MPI_Finalize();
	return 0;
}
