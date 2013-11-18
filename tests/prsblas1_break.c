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
	float tic, toc, toc1;
	float FLOPs;
	float* result;
	float* elapsed;
	float* blas_result;
	float* blas_elapsed;
	int tind;
	int fold;
	int strong;
	float  refsum;
	int dtype;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	n     = read_int(argc, argv, "-n", 1024);
	iters = read_int(argc, argv, "-i", 100);
	fold  = read_int(argc, argv, "-f", 3);
	strong = find_option(argc, argv, "--strong");
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
	for (i =0; i < n; i++) y[i] = 1.0;

	if (rank == 0)
		printf("P N SDOT_RED SDOT_ALL SASUM_RED SASUM_ALL SNRM2 DNR2M_ALL DSUM DSUM_ALL BLAS_SDOT BLAS_SASUM BLAS_SNRM\n");

	// EFFECTIVE FLOP COUNTS: N
	FLOPs = n;
	float MHz = iters * (FLOPs / 1e6);
	
	void *tsum, *gsum;
	F_CARRY_T* ptime;

	tsum = malloc(sizeof(int) + (1 + fold) * sizeof(float) + (1+fold) *sizeof(F_CARRY_T));
	gsum = malloc(sizeof(int) + (1 + fold) * sizeof(float) + (1+fold) *sizeof(F_CARRY_T));

	elapsed = calloc(16 , sizeof(float));
	float* comm = calloc(16 , sizeof(float));
	float* blas_comm = calloc(16 , sizeof(float));
	blas_elapsed = calloc(16 , sizeof(float));
	result = calloc(16 , sizeof(float));
	blas_result = calloc(16 , sizeof(float));

#ifdef CALL_SDOT
	int ione = 1;
	float lsum;
	MPI_Op mpi_nrm2;
	MPI_Op_create(MPI_NRM2, 1, &mpi_nrm2);

	lsum = CALL_SASUM(n, v, ione);
	lsum = CALL_SNRM2(n, v, ione);
	lsum = CALL_SDOT(n, v, ione, y, ione);
#endif

	prsdotI2(MPI_COMM_WORLD, 0, n, v, 1, y, 1, 3, 0, gsum, tsum);
	prsasumI2(MPI_COMM_WORLD, 0, n, v, 1, 3, 0, gsum, tsum);
	prssumI2(MPI_COMM_WORLD, 0, n, v, 1, 3, 0, gsum, tsum);
	prsnrm2I2(MPI_COMM_WORLD, 0, n, v, 1, 3, 0, gsum, tsum);

	for (i = 0; i < iters; i++) {
		// REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		prsdotI2(MPI_COMM_WORLD, 0, n, v, 1, y, 1, 3, 0, gsum, tsum);
		toc = read_timer();
			
		ptime = (F_CARRY_T*) ((char*)gsum + sizeof(int) + fold * (sizeof(float) + sizeof(F_CARRY_T)));
		elapsed[0] += (toc - tic);
		comm[0] += ptime[0];


#ifdef CALL_SDOT
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// local sum
		lsum = CALL_SDOT(n, v, ione, y, ione);

		//REDUCTION
		toc1 = read_timer();
		MPI_Reduce(&lsum, &sum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[0] += (toc - tic);
		blas_comm[0] += (toc - toc1);
		blas_result[0]  = sum;
#endif

//=====
		// REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		prsasumI2(MPI_COMM_WORLD, 0, n, v, 1, 3, 0, gsum, tsum);
		toc = read_timer();

		ptime = (F_CARRY_T*) ((char*)gsum + sizeof(int) + fold * (sizeof(float) + sizeof(F_CARRY_T)));
		elapsed[1] += (toc - tic);
		comm[1] += ptime[0];
			
#ifdef CALL_SDOT
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// local sum
		lsum = CALL_SASUM(n, v, ione);

		//REDUCTION
		toc1 = read_timer();
		MPI_Reduce(&lsum, &sum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[1] += (toc - tic);
		blas_comm[1] += (toc - toc1);
		blas_result[1]  = sum;
#endif

//=====

		// REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		prsnrm2I2(MPI_COMM_WORLD, 0, n, v, 1, 3, 0, gsum, tsum);
		toc = read_timer();
			
		ptime = (F_CARRY_T*) ((char*)gsum + sizeof(int) + (1+fold) * sizeof(float) + fold * sizeof(F_CARRY_T));
		elapsed[2] += (toc - tic);
		comm[2] += ptime[0];

#ifdef CALL_SDOT
		MPI_Barrier(MPI_COMM_WORLD);
		lsum = CALL_SNRM2(n, v, ione);

		tic = read_timer();
		// local sum
		lsum = CALL_SNRM2(n, v, ione);

		//REDUCTION
		toc1 = read_timer();
		MPI_Reduce(&lsum, &sum, 1, MPI_FLOAT, mpi_nrm2, 0, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[2] += (toc - tic);
		blas_comm[2] += (toc - toc1);
		blas_result[2]  = sum;
#endif

//====

		MPI_Barrier(MPI_COMM_WORLD);
		tic = read_timer();
		prssumI2(MPI_COMM_WORLD, 0, n, v, 1, 3, 0, gsum, tsum);
		toc = read_timer();
			
		ptime = (F_CARRY_T*) ((char*)gsum + sizeof(int) + fold * (sizeof(float) + sizeof(F_CARRY_T)));
		elapsed[3] += (toc - tic);
		comm[3] += ptime[0];
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// PRINT OUT
	if (rank == 0) {
		printf("%5d", nprocs);
		printf("%12d", N);
		for (i = 0; i < 4; i++) {
			printf("%8.2f %12.3g %12.3g " , MHz / elapsed[i], elapsed[i], comm[i]);
		}
		for (i = 0; i < 3; i++) {
			printf("%8.2f %12.3g %12.3g" , MHz / blas_elapsed[i], blas_elapsed[i], blas_comm[i]);
		}
		printf("\n");
	}
	
#ifdef CALL_SDOT
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
	free(tsum);
	free(gsum);
	MPI_Finalize();
	return 0;
}
