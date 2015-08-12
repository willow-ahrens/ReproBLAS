/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include "MPI_idxd.h"
#include "../MPI_indexed/MPI_didxd.h"
#include "../types.h"
#include "MPI_reproBLAS.h"
#include "idxdBLAS.h"

void przsumI(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx,
	I_double_Complex* sum
) {
	I_double_Complex local_sum;

	int i;

	RMPI_Init();

	// PEFORM THE LOCAL COMPUTATION
	local_sum = zsumI(N, x, incx);

	zISet(*sum, local_sum);

	// REDUCE THE RESULT TO THE ROOT PROCESSOR
	if (root >= 0)
		MPI_Reduce(&local_sum, sum, 1, MPI_IDOUBLE_COMPLEX,
			MPI_RSUM, root, comm);
	else
		MPI_Allreduce(&local_sum, sum, 1, MPI_IDOUBLE_COMPLEX,
			MPI_RSUM, comm);
}

double complex przsum(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx
) {
	I_double_Complex sum;

	przsumI(comm, root, N, x, incx, &sum);

	return Iconv2z(sum);
}

void przsumI2(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx,
	int fold, double complex* sum,
	double complex* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
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
		myType = MPI_IDOUBLE_COMPLEX;
   	
	// PEFORM THE LOCAL SUM WITH BLOCK SIZE OF 1024
	zsumI1(N, x, incx, fold, local_sum, local_sum + fold);

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

double complex przsum2(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx,
	int fold
) {
	double complex *sum = NULL;
	double complex ret;
	int    me;

	sum = (double complex*) malloc(zISize(fold));

	przsumI2(comm, root, N, x, incx, fold, sum, NULL);

	MPI_Comm_rank(comm, &me);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
	ret = Iconv2z1(fold, sum, sum + fold, 1);

	free(sum);
	return ret;
}

