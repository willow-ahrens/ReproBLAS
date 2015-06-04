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
	float* pin = (float*) invec;
	float* pout = (float*) inoutvec;
	for (i = 0; i < *len; i++) {
		pout[i] = sqrt(pout[i] * pout[i] + pin[i] * pin[i]);
	}
}

int main( int argc, char **argv ) {
	int rank, nprocs;
	int n;
	float complex* v;
	float complex* y;
	int incv;
	float complex sum;
	int i;
	int iters;
	float tic, toc;
	float FLOPs;
	float complex* result;
	float* elapsed;
	float complex* blas_result;
	float* blas_elapsed;
	int tind;
	int fold;
	int check, strong;
	int dtype;
	float K;

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
	K     = read_float(argc, argv, "-K", 1e4);
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
	v     = (float complex*)malloc(n * sizeof(float complex));
	y     = (float complex*)malloc(n * sizeof(float complex));
	v[0]  = K;
	sgenvec(2*n, (float*)v, dtype, 1.0);
	for (i = 0; i < n; i++) y[i] = 1.0;

	elapsed = calloc(16 , sizeof(float));
	blas_elapsed = calloc(16 , sizeof(float));
	result = calloc(16 , sizeof(float complex));
	blas_result = calloc(16 , sizeof(float complex));

	if (rank == 0) {
		printf("%5s", "P");
		printf("%12s", "N");
		printf("%12s", "PRCDOTU_R");
		printf("%12s", "PRCDOTU_A");
		printf("%12s", "PRCDOTC_R");
		printf("%12s", "PRCDOTC_A");
		printf("%12s", "PRSCASUM_R");
		printf("%12s", "PRSCASUM_A");
		printf("%12s", "PRCSUM_R");
		printf("%12s", "PRCSUM_A");
		printf("%12s", "PRSCNRM2_R");
		printf("%12s", "PRSCNRM2_A");
		#ifdef CALL_CDOTU
		printf("%12s", "CDOTU_R");
		printf("%12s", "CDOTU_A");
		#endif
		#ifdef CALL_CDOTC
		printf("%12s", "CDOTC_R");
		printf("%12s", "CDOTC_A");
		#endif
		#ifdef CALL_SCASUM
		printf("%12s", "SCASUM_R");
		printf("%12s", "SCASUM_A");
		#endif
		#ifdef CALL_DNRM2
		printf("%12s", "SCNRM2_R");
		printf("%12s", "SCNRM2_A");
		#endif
	}

	// EFFECTIVE FLOP COUNTS: N
	FLOPs = n;
	float MHz = iters * (FLOPs / 1e6);
	
	float* tsum;

	tsum = (float*) malloc((1 + 2 * fold) * sizeof(float));


	int ione = 1;
	float complex lsum;
	float sumd;
	float gsum;
#ifdef CALL_CDOTU
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
		// CDOTU
		// REDUCE
		ind = 0;
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prcdotu(MPI_COMM_WORLD, 0, n, v, 1, y, 1);
		toc = read_timer();
		
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;

		// ALL REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prcdotu(MPI_COMM_WORLD, -1, n, v, 1, y, 1);
		toc = read_timer();
		
		ind++;
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;

		// CDOTC
		// REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prcdotc(MPI_COMM_WORLD, 0, n, v, 1, y, 1);
		toc = read_timer();
			
		ind++;
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;

		// ALL REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prcdotc(MPI_COMM_WORLD, -1, n, v, 1, y, 1);
		toc = read_timer();
			
		ind++;
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;
#ifdef CALL_CDOTU
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// local sum
		lsum = CALL_CDOTU(n, v, ione, y, ione);

		//REDUCTION
		MPI_Reduce(&lsum, &sum, 1, MPI_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[0] += (toc - tic);
		blas_result [0]  = sum;
#endif

//=====
		// REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prscasum(MPI_COMM_WORLD, 0, n, v, 1);
		toc = read_timer();
			
		ind++;
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;

		// ALL REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prscasum(MPI_COMM_WORLD, -1, n, v, 1);
		toc = read_timer();
			
		ind++;
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;

#ifdef CALL_SCASUM
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// local sum
		sumd = CALL_SCASUM(n, v, ione);

		//REDUCTION
		MPI_Reduce(&sumd, &gsum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		toc  = read_timer();
			
		blas_elapsed[1] += (toc - tic);
		blas_result [1]  = gsum;
#endif

		// REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prcsum(MPI_COMM_WORLD, 0, n, v, 1);
		toc = read_timer();
			
		ind++;
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;

		// ALL REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prcsum(MPI_COMM_WORLD, -1, n, v, 1);
		toc = read_timer();
			
		ind++;
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;
//=====
		// REDUCE
		MPI_Barrier(MPI_COMM_WORLD);
		tic = read_timer();
		sum = prscnrm2(MPI_COMM_WORLD, 0, n, v, 1);
		toc = read_timer();
			
		ind++;
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;

		// ALL REDUCE
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		sum = prscnrm2(MPI_COMM_WORLD, -1, n, v, 1);
		toc = read_timer();
			
		ind++;
		elapsed[ind] += (toc - tic);
		result [ind]  = sum;

#ifdef CALL_SCNRM2
		MPI_Barrier(MPI_COMM_WORLD);

		tic = read_timer();
		// local sum
		sumd = CALL_DNRM2(n, v, ione);

		//REDUCTION
		MPI_Reduce(&sumd, &gsum, 1, MPI_FLOAT, mpi_nrm2, 0, MPI_COMM_WORLD);
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
		float complex* lv = (float complex*) malloc((int)N * sizeof(float complex));
		float complex* ly = (float complex*) malloc((int)N * sizeof(float complex));
		float err, rerr;
		if (rank == 0)
			printf("%17s", "Reproducibility");

		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Allgather(v, n, MPI_COMPLEX, lv, n, MPI_COMPLEX, MPI_COMM_WORLD);
		MPI_Allgather(y, n, MPI_COMPLEX, ly, n, MPI_COMPLEX, MPI_COMM_WORLD);

//#ifdef CALL_PRCDOT
		ind = 0;
		sum = rcdotu(N, lv, 1, ly, 1);
		err = RELERR(creal(sum), creal(result[ind+1]));
		MPI_Reduce(&err, &rerr, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			printf("%12.2g", RELERR(creal(sum), creal(result[ind])));
			printf("%12.2g", rerr);
		}		

		ind += 2;
		sum = rcdotc(N, lv, 1, ly, 1);
		err = RELERR(creal(sum), creal(result[ind + 1]));
		MPI_Reduce(&err, &rerr, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			printf("%12.2g", RELERR(creal(sum), creal(result[ind])));
			printf("%12.2g", rerr);
		}		
//#endif

//#ifdef CALL_PRCASUM
		sumd = rscasum(N, lv, 1);
		ind += 2;
		err = RELERR(sumd, creal(result[ind + 1]));
		MPI_Reduce(&err, &rerr, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			printf("%12.2g", RELERR(sumd, creal(result[ind])));
			printf("%12.2g", rerr);
		}		
//#endif

//#ifdef CALL_PRCSUM
		sum = rcsum(N, lv, 1);

		ind += 2;
		err = RELERR(creal(sum), creal(result[ind+1]));
		MPI_Reduce(&err, &rerr, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			printf("%12.2g", RELERR(creal(sum), creal(result[ind])));
			printf("%12.2g", rerr);
		}		
//#endif

//#ifdef CALL_PRCNRM2
		sumd = rscnrm2(N, lv, 1);
		ind += 2;
		err = RELERR(sumd, creal(result[ind+1]));
		MPI_Reduce(&err, &rerr, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			printf("%12.2g", RELERR(sumd, creal(result[ind])));
			printf("%12.2g", rerr);
		}		
//#endif

//=====
#ifdef CALL_CDOTU

		sum = CALL_CDOTU(N, lv, ione, ly, ione);
		ind += 2;
		err = RELERR(creal(sum), creal(blas_result[1]));
		MPI_Reduce(&err, &rerr, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

		if (rank == 0) {
			printf("%12.2g", RELERR(creal(sum), creal(blas_result[0])));
			printf("%12.2g", rerr);
		}
#endif

#ifdef CALL_SCASUM
		sumd = CALL_SCASUM(N, lv, ione);
		err = RELERR(sum, creal(blas_result[3]));
		MPI_Reduce(&err, &rerr, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

		if (rank == 0) {
			printf("%12.2g", RELERR(creal(sumd, creal(blas_result[2])));
			printf("%12.2g", rerr);
		}
#endif

#ifdef CALL_DNRM2
		sum = CALL_SCNRM2(N, lv, ione);

		err = RELERR(creal(sum), creal(blas_result[5]));
		MPI_Reduce(&err, &rerr, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		
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
