/*
 *  Created   13/10/25   H.D. Nguyen
 */

//#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include "MPI_indexed.h"
#include "../MPI_indexed/MPI_dindexed.h"
#include "MPI_reproBLAS.h"
#include "indexedBLAS.h"

void prddotI(
	MPI_Comm comm, int root,
	int N,
	double* x, int incx,
	double* y, int incy,
	I_double* sum
) {
	Idouble local_sum;
	int i;

	RMPI_Init();
	// PEFORM THE LOCAL SUM WITH BLOCK SIZE OF 1024
	dISetZero(local_sum);
	local_sum = ddotI(N, x, incx, y, incy);

	dISet(*sum, local_sum);
	// REDUCE THE RESULT TO THE ROOT PROCESSOR
	if (root >= 0)
		MPI_Reduce(&local_sum, sum, 1, MPI_IDOUBLE, MPI_RSUM, root, comm);
	else
		MPI_Allreduce(&local_sum, sum, 1, MPI_IDOUBLE, MPI_RSUM, comm);
}

double prddot(
	MPI_Comm comm, int root,
	int N,
	double* x, int incx,
	double* y, int incy
) {
	Idouble sum;
	int i;

	prddotI(comm, root, N, x, incx, y, incy, &sum);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
	return Iconv2d(sum);
}

void prddotI2(
	MPI_Comm comm, int root,
	int N,
	double* x, int incx,
	double* y, int incy,
	int fold, int W, double* sum,
	double* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
) {
	int me, i;

	MPI_Datatype myType;
	int created = 0;
	if (local_sum == NULL) {
		local_sum = (double*) malloc((2 * fold) * sizeof(double));
		created = 1;
	}
	for (i =0; i < 2 * fold; i++) local_sum[i] = 0.0;

	// REGISTER THE MPI TYPE FOR THE EXTENDED INDEXED FP NUMBER
	if (fold != DEFAULT_FOLD)
		dIMPICreate(fold, &myType);
	else
		myType = MPI_IDOUBLE;
   	
	// PEFORM THE LOCAL SUM WITH BLOCK SIZE OF 1024
	ddotI1(N, x, incx, y, incy, fold, W, local_sum, local_sum + fold);

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

double prddot2(
	MPI_Comm comm, int root,
	int N,
	double* x, int incx,
	double* y, int incy,
	int fold, int W 
) {
	double *sum = NULL;
	double ret;
	int    me;

	sum = (double*) malloc((2 * fold) * sizeof(double));

	prddotI2(comm, root, N, x, incx, y, incy, fold, W, sum, NULL);

	MPI_Comm_rank(comm, &me);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
	ret = Iconv2d1(fold, sum, sum + fold, 1);

	free(sum);
	return ret;
}

