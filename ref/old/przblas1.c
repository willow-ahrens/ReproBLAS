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
	double complex* v;
	double complex* y;
	int incv;
	double complex sum;
	int i;
	int iters;
	double tic, toc;
	double FLOPs;
	double complex* result;
	double* elapsed;
	double complex* blas_result;
	double* blas_elapsed;
	int tind;
	int fold;
	int check, strong;
	int dtype;
	double K;

	MPI_Init(&argc, &argv);
	RMPI_Init();
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
	v     = (double complex*)malloc(n * sizeof(double complex));
	y     = (double complex*)malloc(n * sizeof(double complex));
	v[0]  = K;
	dgenvec(2*n, (double*)v, dtype, 1.0);
	for (i = 0; i < n; i++) y[i] = 1.0;

	elapsed = calloc(16 , sizeof(double));
	blas_elapsed = calloc(16 , sizeof(double));
	result = calloc(16 , sizeof(double complex));
	blas_result = calloc(16 , sizeof(double complex));

	if (rank == 0) {
		printf("%5s", "P");
		printf("%12s", "N");
		printf("%12s", "PRZDOTU_R");
		printf("%12s", "PRZDOTU_A");
		printf("%12s", "PRZDOTC_R");
		printf("%12s", "PRZDOTC_A");
		printf("%12s", "PRDZASUM_R");
		printf("%12s", "PRDZASUM_A");
		printf("%12s", "PRZSUM_R");
		printf("%12s", "PRZSUM_A");
		printf("%12s", "PRDZNRM2_R");
		printf("%12s", "PRDZNRM2_A");
		#ifdef CALL_ZDOTU
		printf("%12s", "ZDOTU_R");
		printf("%12s", "ZDOTU_A");
		#endif
		#ifdef CALL_ZDOTC
		printf("%12s", "ZDOTC_R");
		printf("%12s", "ZDOTC_A");
		#endif
		#ifdef CALL_DZASUM
		printf("%12s", "DZASUM_R");
		printf("%12s", "DZASUM_A");
		#endif
		#ifdef CALL_DNRM2
		printf("%12s", "DZNRM2_R");
		printf("%12s", "DZNRM2_A");
		#endif
	}

	// EFFECTIVE FLOP COUNTS: N
	FLOPs = n;
	double MHz = iters * (FLOPs / 1e6);
	
	double* tsum;

	tsum = (double*) malloc((1 + 2 * fold) * sizeof(double));


	int ione = 1;
	double complex lsum;
	double sumd;
	double gsum;
#ifdef CALL_ZDOTU
	MPI_Op mpi_nrm2;
	MPI_Op mpi_sum;
	MPI_Op_create(MPI_NRM2, 1, &mpi_nrm2);
	if (rank == 0)
	printf(" BLAS_DDOT BLAS_DASUM BLAS_DNRM");
#endif
	if (rank == 0)
	printf("\n");
	int ind;

	for (i = 0; i < iters; i++) {
		// ZDOTU
		// REDUCE
		ind = 0;
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = przdotu(MPI_COMM_WORLD, 0, n, v, 1, y, 1);
		toc = read_timer();
		
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;

		// ALL REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = przdotu(MPI_COMM_WORLD, -1, n, v, 1, y, 1);
		toc = read_timer();
		
		ind++;
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;

		// ZDOTC
		// REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = przdotc(MPI_COMM_WORLD, 0, n, v, 1, y, 1);
		toc = read_timer();
			
		ind++;
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;

		// ALL REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = przdotc(MPI_COMM_WORLD, -1, n, v, 1, y, 1);
		toc = read_timer();
			
		ind++;
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;
#ifdef CALL_ZDOTU
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// local sum
		lsum = CALL_ZDOTU(n, v, ione, y, ione);

		//REDUCTION
		MPI_Reduce(&lsum, &sum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[0] += (toc - tic);
		blas_result [0]  = sum;
#endif

//=====
		// REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prdzasum(MPI_COMM_WORLD, 0, n, v, 1);
		toc = read_timer();
			
		ind++;
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;

		// ALL REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prdzasum(MPI_COMM_WORLD, -1, n, v, 1);
		toc = read_timer();
			
		ind++;
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;

#ifdef CALL_DZASUM
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// local sum
		sumd = CALL_DZASUM(n, v, ione);

		//REDUCTION
		MPI_Reduce(&sumd, &gsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[1] += (toc - tic);
		blas_result [1]  = gsum;
#endif

		// REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = przsum(MPI_COMM_WORLD, 0, n, v, 1);
		toc = read_timer();
			
		ind++;
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;

		// ALL REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = przsum(MPI_COMM_WORLD, -1, n, v, 1);
		toc = read_timer();
			
		ind++;
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;
//=====
		// REDUCE
		MPI_Barrier(MPI_COMM_WORLD);
		tic = read_timer();
		sum = prdznrm2(MPI_COMM_WORLD, 0, n, v, 1);
		toc = read_timer();
			
		ind++;
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;

		// ALL REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prdznrm2(MPI_COMM_WORLD, -1, n, v, 1);
		toc = read_timer();
			
		ind++;
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;

#ifdef CALL_DZNRM2
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// local sum
		sumd = CALL_DNRM2(n, v, ione);

		//REDUCTION
		MPI_Reduce(&sumd, &gsum, 1, MPI_DOUBLE, mpi_nrm2, 0, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[2] += (toc - tic);
		blas_result [2]  = sum;
#endif

//====

	}

	MPI_Barrier(MPI_COMM_WORLD);

	// PRINT OUT
	if (rank == 0) {
		printf("%5d", nprocs);
		printf("%12d", N);
		for (i = 0; i < 10; i++) {
			printf("%12.2g" , MHz / elapsed[i]);
		}
#ifdef CALL_DDOT
		for (i = 0; i < 3; i++) {
			printf("%12.2g" , MHz / blas_elapsed[i]);
		}
#endif
		printf("\n");
	}

	if (check) {
		double complex* lv = (double complex*) malloc((int)N * sizeof(double complex));
		double complex* ly = (double complex*) malloc((int)N * sizeof(double complex));
		double err, rerr;
		if (rank == 0)
			printf("%17s", "Reproducibility");

		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Allgather(v, n, MPI_DOUBLE_COMPLEX, lv, n, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
		MPI_Allgather(y, n, MPI_DOUBLE_COMPLEX, ly, n, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);

//#ifdef CALL_PRZDOT
		ind = 0;
		sum = rzdotu(N, lv, 1, ly, 1);
		err = RELERR(creal(sum), creal(result[ind+1]));
		MPI_Reduce(&err, &rerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			printf("%12.2g", RELERR(creal(sum), creal(result[ind])));
			printf("%12.2g", rerr);
		}		

		ind += 2;
		sum = rzdotc(N, lv, 1, ly, 1);
		err = RELERR(creal(sum), creal(result[ind + 1]));
		MPI_Reduce(&err, &rerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			printf("%12.2g", RELERR(creal(sum), creal(result[ind])));
			printf("%12.2g", rerr);
		}		
//#endif

//#ifdef CALL_PRZASUM
		sumd = reproBLAS_rdzasum(N, lv, 1);
		ind += 2;
		err = RELERR(sumd, creal(result[ind + 1]));
		MPI_Reduce(&err, &rerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			printf("%12.2g", RELERR(sumd, creal(result[ind])));
			printf("%12.2g", rerr);
		}		
//#endif

//#ifdef CALL_PRZSUM
		sum = rzsum(N, lv, 1);

		ind += 2;
		err = RELERR(creal(sum), creal(result[ind+1]));
		MPI_Reduce(&err, &rerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			printf("%12.2g", RELERR(creal(sum), creal(result[ind])));
			printf("%12.2g", rerr);
		}		
//#endif

//#ifdef CALL_PRZNRM2
		sumd = reproBLAS_rdznrm2(N, lv, 1);
		ind += 2;
		err = RELERR(sumd, creal(result[ind+1]));
		MPI_Reduce(&err, &rerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			printf("%12.2g", RELERR(sumd, creal(result[ind])));
			printf("%12.2g", rerr);
		}		
//#endif

//=====
#ifdef CALL_ZDOTU

		sum = CALL_ZDOTU(N, lv, ione, ly, ione);
		ind += 2;
		err = RELERR(creal(sum), creal(blas_result[1]));
		MPI_Reduce(&err, &rerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		if (rank == 0) {
			printf("%12.2g", RELERR(creal(sum), creal(blas_result[0])));
			printf("%12.2g", rerr);
		}
#endif

#ifdef CALL_DZASUM
		sumd = CALL_DZASUM(N, lv, ione);
		err = RELERR(sum, creal(blas_result[3]));
		MPI_Reduce(&err, &rerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		if (rank == 0) {
			printf("%12.2g", RELERR(creal(sumd, creal(blas_result[2])));
			printf("%12.2g", rerr);
		}
#endif

#ifdef CALL_DNRM2
		sum = CALL_DZNRM2(N, lv, ione);

		err = RELERR(creal(sum), creal(blas_result[5]));
		MPI_Reduce(&err, &rerr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		if (rank == 0) {
			printf("%12.2g", RELERR(creal(sum), creal(blas_result[4])));
			printf("%12.2g", rerr);
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
