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
	float* pin = (float*) invec;
	float* pout = (float*) inoutvec;
	for (i = 0; i < *len; i++) {
		pout[i] = sqrt(pout[i] * pout[i] + pin[i] * pin[i]);
	}
}

int main( int argc, char **argv ) {
	int rank, nprocs;
	int n;
	float* v;
	float* y;
	int incv;
	float sum;
	int i;
	int iters;
	float tic, toc;
	float FLOPs;
	float* result;
	float* elapsed;
	float* blas_result;
	float* blas_elapsed;
	int tind;
	int fold;
	int strong;
	int check;
	float  refsum;
	int dtype;
	float err, rerr;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	n     = read_int(argc, argv, "-n", 1024);
	iters = read_int(argc, argv, "-i", 100);
	fold  = read_int(argc, argv, "-f", 3);
	strong = find_option(argc, argv, "--strong");
	check = find_option(argc, argv, "--check") >= 0;
	dtype = read_int(argc, argv, "-d", 1);
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
	v     = (float*)malloc(n * sizeof(float));
	y     = (float*)malloc(n * sizeof(float));
	sgenvec(n, v, dtype, 1.0);
	for (i = 0; i < n; i++) y[i] = 1.0;

	if (rank == 0) {
		printf("%5s", "P");
		printf("%12s", "N");
		#ifdef CALL_PRSDOT
		printf("%10s", "PRSDOT_R");
		printf("%10s", "PRSDOT_A");
		#endif
		#ifdef CALL_PRSASUM
		printf("%10s", "PRSASUM_R");
		printf("%10s", "PRSASUM_A");
		#endif
		#ifdef CALL_PRSSUM
		printf("%10s", "PRSSUM_R");
		printf("%10s", "PRSSUM_A");
		#endif
		#ifdef CALL_PRSNRM2
		printf("%10s", "PRSNRM2_R");
		printf("%10s", "PRSNRM2_A");
		#endif
		#ifdef CALL_SDOT
		printf("%10s", "SDOT_R");
		printf("%10s", "SDOT_A");
		#endif
		#ifdef CALL_SASUM
		printf("%10s", "SASUM_R");
		printf("%10s", "SASUM_A");
		#endif
		#ifdef CALL_SNRM2
		printf("%10s", "SNRM2_R");
		printf("%10s", "SNRM2_A");
		#endif
		printf("\n");
	}

	// EFFECTIVE FLOP COUNTS: N
	FLOPs = n;
	float MHz = iters * (FLOPs / 1e6);
	
	float* tsum;

	tsum = (float*) malloc((1 + 2 * fold) * sizeof(float));

	elapsed = calloc(16 , sizeof(float));
	blas_elapsed = calloc(16 , sizeof(float));
	result = calloc(16 , sizeof(float));
	blas_result = calloc(16 , sizeof(float));

	int ione = 1;
	float lsum;

#ifdef CALL_SNRM2
	MPI_Op mpi_nrm2;
	MPI_Op_create(MPI_NRM2, 1, &mpi_nrm2);
#endif

	for (i = 0; i < iters; i++) {
		// REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

#ifdef CALL_PRSDOT
		tic = read_timer();
		sum = prsdot(MPI_COMM_WORLD, 0, n, v, 1, y, 1);
		toc = read_timer();
			
		elapsed[0] += (toc - tic);
		result [0]  = sum;

		// ALL REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prsdot(MPI_COMM_WORLD, -1, n, v, 1, y, 1);
		toc = read_timer();
			
		elapsed[1] += (toc - tic);
		result [1]  = sum;

#endif

#ifdef CALL_SDOT
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// local sum
		lsum = CALL_SDOT(n, v, ione, y, ione);

		//REDUCTION
		MPI_Reduce(&lsum, &sum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[0] += (toc - tic);
		blas_result [0]  = sum;

		//--- ALL REUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// local sum
		lsum = CALL_SDOT(n, v, ione, y, ione);

		//REDUCTION
		MPI_Allreduce(&lsum, &sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[1] += (toc - tic);
		blas_result [1]  = sum;

#endif

//=====
		// REDUCE
#ifdef CALL_PRSASUM
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prsasum(MPI_COMM_WORLD, 0, n, v, 1);
		toc = read_timer();
			
		elapsed[2] += (toc - tic);
		result [2]  = sum;

		// ALL REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prsasum(MPI_COMM_WORLD, -1, n, v, 1);
		toc = read_timer();
			
		elapsed[3] += (toc - tic);
		result [3]  = sum;
#endif

#ifdef CALL_PRSSUM
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prssum(MPI_COMM_WORLD, 0, n, v, 1);
		toc = read_timer();
			
		elapsed[6] += (toc - tic);
		result [6]  = sum;

		// ALL REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prssum(MPI_COMM_WORLD, -1, n, v, 1);
		toc = read_timer();
			
		elapsed[7] += (toc - tic);
		result [7]  = sum;
#endif

#ifdef CALL_SASUM
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// local sum
		lsum = CALL_SASUM(n, v, ione);

		//REDUCTION
		MPI_Reduce(&lsum, &sum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[2] += (toc - tic);
		blas_result [2]  = sum;
		//--- all reduce ---

		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// local sum
		lsum = CALL_SASUM(n, v, ione);

		//REDUCTION
		MPI_Allreduce(&lsum, &sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[3] += (toc - tic);
		blas_result [3]  = sum;
#endif

//=====
#ifdef CALL_PRSNRM2
		// REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prsnrm2(MPI_COMM_WORLD, 0, n, v, 1);
		toc = read_timer();
			
		elapsed[4] += (toc - tic);
		result [4]  = sum;

		// ALL REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prsnrm2(MPI_COMM_WORLD, -1, n, v, 1);
		toc = read_timer();
			
		elapsed[5] += (toc - tic);
		result [5]  = sum;
#endif

#ifdef CALL_SNRM2
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// local sum
		lsum = CALL_SNRM2(n, v, ione);

		//REDUCTION
		MPI_Reduce(&lsum, &sum, 1, MPI_FLOAT, mpi_nrm2, 0, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[4] += (toc - tic);
		blas_result [4]  = sum;
		//---

		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// local sum
		lsum = CALL_SNRM2(n, v, ione);

		//REDUCTION
		MPI_Allreduce(&lsum, &sum, 1, MPI_FLOAT, mpi_nrm2, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[5] += (toc - tic);
		blas_result [5]  = sum;
#endif
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// PRINT OUT
	if (rank == 0) {
		printf("%5d", nprocs);
		printf("%12d", N);
		/*
		for (i = 0; i < 8; i++) {
			printf("%12.2g" , MHz / elapsed[i]);
		}
		for (i = 0; i < 3; i++) {
			printf("%12.2g" , MHz / blas_elapsed[i]);
		}
		*/

		#ifdef CALL_PRSDOT
		printf("%10.2g" , MHz / elapsed[0]);
		printf("%10.2g" , MHz / elapsed[1]);
		#endif

		#ifdef CALL_PRSASUM
		printf("%10.2g" , MHz / elapsed[2]);
		printf("%10.2g" , MHz / elapsed[3]);
		#endif

		#ifdef CALL_PRSSUM
		printf("%10.2g" , MHz / elapsed[2]);
		printf("%10.2g" , MHz / elapsed[3]);
		#endif

		#ifdef CALL_PRSNRM2
		printf("%10.2g" , MHz / elapsed[4]);
		printf("%10.2g" , MHz / elapsed[5]);
		#endif

		#ifdef CALL_SDOT
		printf("%10.2g" , MHz / blas_elapsed[0]);
		printf("%10.2g" , MHz / blas_elapsed[1]);
		#endif

		#ifdef CALL_SASUM
		printf("%10.2g" , MHz / blas_elapsed[2]);
		printf("%10.2g" , MHz / blas_elapsed[3]);
		#endif

		#ifdef CALL_SNRM2
		printf("%10.2g" , MHz / blas_elapsed[4]);
		printf("%10.2g" , MHz / blas_elapsed[5]);
		#endif
		printf("\n");
	}


	if (check) {
		float* lv = (float*) malloc((int)N * sizeof(float));
		float* ly = (float*) malloc((int)N * sizeof(float));
		if (rank == 0)
			printf("%17s", "Reproducibility");

		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Allgather(v, n, MPI_FLOAT, lv, n, MPI_FLOAT, MPI_COMM_WORLD);
		MPI_Allgather(y, n, MPI_FLOAT, ly, n, MPI_FLOAT, MPI_COMM_WORLD);

#ifdef CALL_PRSDOT
		sum = rsdot(N, lv, 1, ly, 1);
		err = RELERR(sum, result[1]);
		MPI_Reduce(&err, &rerr, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			printf("%10.2g", RELERR(sum, result[0]));
			printf("%10.2g", rerr);
		}		
#endif

#ifdef CALL_PRSASUM
		sum = rsasum(N, lv, 1);

		err = RELERR(sum, result[3]);
		MPI_Reduce(&err, &rerr, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			printf("%10.2g", RELERR(sum, result[2]));
			printf("%10.2g", rerr);
		}		
#endif

#ifdef CALL_PRSSUM
		sum = rssum(N, lv, 1);

		err = RELERR(sum, result[7]);
		MPI_Reduce(&err, &rerr, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			printf("%10.2g", RELERR(sum, result[6]));
			printf("%10.2g", rerr);
		}		
#endif

#ifdef CALL_PRSNRM2
		sum = rsnrm2(N, lv, 1);
		err = RELERR(sum, result[5]);
		MPI_Reduce(&err, &rerr, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			printf("%10.2g", RELERR(sum, result[4]));
			printf("%10.2g", rerr);
		}		
#endif

//=====
#ifdef CALL_SDOT

		sum = CALL_SDOT(N, lv, ione, ly, ione);
		err = RELERR(sum, blas_result[1]);
		MPI_Reduce(&err, &rerr, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

		if (rank == 0) {
			printf("%10.2g", RELERR(sum, blas_result[0]));
			printf("%10.2g", rerr);
		}
#endif

#ifdef CALL_SASUM
		sum = CALL_SASUM(N, lv, ione);
		err = RELERR(sum, blas_result[3]);
		MPI_Reduce(&err, &rerr, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

		if (rank == 0) {
			printf("%10.2g", RELERR(sum, blas_result[2]));
			printf("%10.2g", rerr);
		}
#endif

#ifdef CALL_SNRM2
		sum = CALL_SNRM2(N, lv, ione);

		err = RELERR(sum, blas_result[5]);
		MPI_Reduce(&err, &rerr, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		
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


#ifdef CALL_SNRM2
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
