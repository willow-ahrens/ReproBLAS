/*
 *  Created   13/10/25   H.D. Nguyen
 */

//#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "MPI_indexed.h"
#include "../MPI_indexed/MPI_dindexed.h"
#include "MPI_reproBLAS.h"
#include "indexedBLAS.h"

void prdznrm2I(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx,
	double* sum
) {
	double local_sum[1 + 2 * DEFAULT_FOLD];
	int i;
	RMPI_Init();

	for (i = 0; i < 1 + 2 * DEFAULT_FOLD; i++)
		local_sum[i] = 0.0;
	local_sum[0] = dznrm2I1(N, x, incx, DEFAULT_FOLD, 0, (double*)local_sum + 1, local_sum + 1 + DEFAULT_FOLD);

	for (i = 0; i < 1 + 2 * DEFAULT_FOLD; i++) {
		sum[i] = local_sum[i];
	}

	if (root >= 0)
		MPI_Reduce(local_sum,sum,1,MPI_IDOUBLE_SCALE, MPI_RNRM2, root, comm);
	else
		MPI_Allreduce(local_sum,sum, 1,MPI_IDOUBLE_SCALE, MPI_RNRM2, comm);
}

double prdznrm2(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx
) {
	double sum[1 + 2 * DEFAULT_FOLD];

	prdznrm2I(comm, root, N, x, incx, sum);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
	return sum[0] * sqrt(Iconv2d1(DEFAULT_FOLD, sum+1, sum + 1 + DEFAULT_FOLD, 1));
}


void prdznrm2I2(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx,
	int fold, int W, double* sum,
	double* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
) {
	int me, i;

	MPI_Datatype myType;

	int created = 0;
	if (local_sum == NULL) {
		local_sum = (double*) malloc(sizeof(double) + dISize(fold));
		created = 1;
	}

	for (i = 0; i < 1 + 2 * fold; i++)
		local_sum[i] = 0.0;

	// REGISTER THE MPI TYPE FOR THE EXTENDED INDEXED FP NUMBER
	if (fold != DEFAULT_FOLD)
		dIMPIScaleCreate(fold, &myType);
	else
		myType = MPI_IDOUBLE_SCALE;
   	
	// PEFORM THE LOCAL SUM WITH BLOCK SIZE OF 1024
	local_sum[0] = dznrm2I1(N, x, incx, fold, W, local_sum + 1, local_sum + 1 + fold);

	// REDUCE THE RESULT TO THE ROOT PROCESSOR
	if (root >= 0)
		MPI_Reduce(local_sum, sum, 1, myType, MPI_RNRM2, root, comm);
	else
		MPI_Allreduce(local_sum, sum, 1, myType, MPI_RNRM2, comm);
			
	// DESTROY THE REGISTERED MPI TYPE
	if (fold != DEFAULT_FOLD)
		MPI_Type_free(&myType);

	if (created) {
		free(local_sum);
	}
}

double prdznrm22(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx,
	int fold, int W
) {
	double *sum = NULL;
	double *local_sum = NULL;
	double ret;
	int    me;

	sum = (double*) malloc(sizeof(double) + dISize(fold));

	prdznrm2I2(comm, root, N, x, incx, fold, W, sum, NULL);

	MPI_Comm_rank(comm, &me);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
	ret = sum[0] * sqrt(Iconv2d1(fold, sum + 1, sum + 1 + fold, 1));

	free(sum);
	return ret;
}

