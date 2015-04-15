/*
 *  Created   13/10/25   H.D. Nguyen
 */

//#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include "MPI_indexed.h"
#include "../MPI_indexed/MPI_dindexed.h"
#include "../types.h"
#include "MPI_reproBLAS.h"
#include "indexedBLAS.h"

void przdotI(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx,
	double complex* y, int incy,
	I_double_Complex* sum,
	int conj
) {
	I_double_Complex local_sum;
	int i;

	RMPI_Init();
	// PEFORM THE LOCAL SUM WITH BLOCK SIZE OF 1024

	if (conj == 1)
		local_sum = zdotcI(N, x, incx, y, incy);
	else
		local_sum = zdotuI(N, x, incx, y, incy);

	zISet(*sum, local_sum);

	// REDUCE THE RESULT TO THE ROOT PROCESSOR
	if (root >= 0)
		MPI_Reduce(&local_sum, sum, 1, MPI_IDOUBLE_COMPLEX,
			MPI_RSUM, root, comm);
	else
		MPI_Allreduce(&local_sum, sum, 1, MPI_IDOUBLE_COMPLEX,
			MPI_RSUM, comm);
}

void przdotcI(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx,
	double complex* y, int incy,
	I_double_Complex* sum
) {
	przdotI(comm, root, N, x, incx, y, incy, sum, 1);
}
void przdotuI(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx,
	double complex* y, int incy,
	I_double_Complex* sum
) {
	przdotI(comm, root, N, x, incx, y, incy, sum, 0);
}

double complex przdotc(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx,
	double complex* y, int incy
) {
	I_double_Complex sum;
	int i;

	przdotI(comm, root, N, x, incx, y, incy, &sum, 1);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
	return Iconv2z(sum);
}

double complex przdotu(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx,
	double complex* y, int incy
) {
	I_double_Complex sum;
	int i;

	przdotI(comm, root, N, x, incx, y, incy, &sum, 0);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
	return Iconv2z(sum);
}


void przdotI2(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx,
	double complex* y, int incy,
	int fold, double complex* sum,
	double complex* local_sum,	// WORKING BUFFER TO STORE LOCAL SUM
	int conj
) {
	int me, i;

	MPI_Datatype myType;
	int created = 0;
	if (local_sum == NULL) {
		local_sum = (double complex*) malloc(zISize(fold));
		created = 1;
	}

	for (i = 0; i < 2 * fold; i++)
		local_sum[i] = 0.0;

	// REGISTER THE MPI TYPE FOR THE EXTENDED INDEXED FP NUMBER
	if (fold != DEFAULT_FOLD)
		zIMPICreate(fold, &myType);
	else
		myType = MPI_IDOUBLE;
   	
	// PEFORM THE LOCAL SUM WITH BLOCK SIZE OF 1024
	if (conj == 1)
		zdotcI1(N, x, incx, y, incy, fold, local_sum, local_sum + fold);
	else
		zdotuI1(N, x, incx, y, incy, fold, local_sum, local_sum + fold);

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

void przdotcI2(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx,
	double complex* y, int incy,
	int fold, double complex* sum,
	double complex* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
) {
	przdotI2(comm, root, N, x, incx, y, incy, fold, sum, local_sum, 1);
}
void przdotuI2(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx,
	double complex* y, int incy,
	int fold, double complex* sum,
	double complex* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
) {
	przdotI2(comm, root, N, x, incx, y, incy, fold, sum, local_sum, 0);
}

double complex przdotc2(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx,
	double complex* y, int incy,
	int fold
) {
	double complex *sum = NULL;
	double complex ret;
	int    me;

	sum = (double complex*) malloc(zISize(fold));

	przdotcI2(comm, root, N, x, incx, y, incy, fold, sum, NULL);

	MPI_Comm_rank(comm, &me);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
	ret = Iconv2z1(fold, sum, sum + fold, 1);

	free(sum);
	return ret;
}

double complex przdotu2(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx,
	double complex* y, int incy,
	int fold
) {
	double complex *sum = NULL;
	double complex ret;
	int    me;

	sum = (double complex*) malloc(zISize(fold));

	przdotuI2(comm, root, N, x, incx, y, incy, fold, sum, NULL);

	MPI_Comm_rank(comm, &me);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
	ret = Iconv2z1(fold, sum, sum + fold, 1);

	free(sum);
	return ret;
}

