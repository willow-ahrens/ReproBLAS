/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include "MPI_indexed.h"
#include "../MPI_indexed/MPI_sindexed.h"
#include "MPI_reproBLAS.h"
#include "indexedBLAS.h"
#include "../types.h"

void prcsumI(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx,
	I_float_Complex* sum
) {
	I_float_Complex local_sum;

	int i;

	RMPI_Init();

	// PEFORM THE LOCAL COMPUTATION
	local_sum = csumI(N, x, incx);

	cISet(*sum, local_sum);

	// REDUCE THE RESULT TO THE ROOT PROCESSOR
	if (root >= 0)
		MPI_Reduce(&local_sum, sum, 1, MPI_ICOMPLEX,
			MPI_RSUM, root, comm);
	else
		MPI_Allreduce(&local_sum, sum, 1, MPI_ICOMPLEX,
			MPI_RSUM, comm);
}

float complex prcsum(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx
) {
	I_float_Complex sum;

	prcsumI(comm, root, N, x, incx, &sum);

	return Iconv2c(sum);
}

void prcsumI2(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx,
	int fold, int W, float complex* sum,
	float complex* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
) {
	int me, i;

	MPI_Datatype myType;
	int created = 0;
	if (local_sum == NULL) {
		local_sum = (float complex*) malloc(cISize(fold));
		created = 1;
	}

	float* C = (float*) (local_sum + fold);
	for (i = 0; i < fold; i++) {
		CSET_(local_sum[i], 0.0, 0.0);
		C[2*i] = C[2*i+1] = 0;
	}

	// REGISTER THE MPI TYPE FOR THE EXTENDED INDEXED FP NUMBER
	if (fold != DEFAULT_FOLD)
		sIMPICreate(2*fold, &myType);
	else
		myType = MPI_ICOMPLEX;
   	
	// PEFORM THE LOCAL SUM WITH BLOCK SIZE OF 1024
	csumI1(N, x, incx, fold, W, local_sum, C);

	// REDUCE THE RESULT TO THE ROOT PROCESSOR
	if (root >= 0)
		MPI_Reduce(local_sum, sum, 1, myType, MPI_RSUM, root, comm);
	else
		MPI_Allreduce(local_sum, sum, 1, myType, MPI_RSUM, comm);
			
	// DESTROY THE REGISTERED MPI TYPE
	if (fold != DEFAULT_FOLD)
		MPI_Type_free(&myType);

	if (created) {
		free(local_sum);
	}
}

float complex prcsum2(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx,
	int fold, int W 
) {
	float complex *sum = NULL;
	float complex ret;
	int    me;

	sum = (float complex*) malloc(cISize(fold));

	prcsumI2(comm, root, N, x, incx, fold, W, sum, NULL);

	MPI_Comm_rank(comm, &me);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
	ret = Iconv2c1(fold, sum, (float*)(sum + fold), 1);

	free(sum);
	return ret;
}

