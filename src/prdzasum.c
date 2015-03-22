/*
 *  Created   13/10/25   H.D. Nguyen
 */

//#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include "rblas1.h"
#include "MPI_Indexed.h"
#include "MPI_Indexed/MPI_dIndexed.h"
#include "prblas.h"

void prdzasumI(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx,
	I_double* sum
) {
	Idouble local_sum;
	int i;

	RMPI_Init();
	// PEFORM THE LOCAL SUM WITH BLOCK SIZE OF 1024
	local_sum = dzasumI(N, x, incx);

	dISet(*sum, local_sum);

	// REDUCE THE RESULT TO THE ROOT PROCESSOR
	if (root >= 0)
		MPI_Reduce(&local_sum, sum, 1, MPI_IDOUBLE, MPI_RSUM, root, comm);
	else
		MPI_Allreduce(&local_sum, sum, 1, MPI_IDOUBLE, MPI_RSUM, comm);
}

double prdzasum(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx
) {
	Idouble sum;
	int    me;

	prdzasumI(comm, root, N, x, incx, &sum);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
	return Iconv2d(sum);
}

void prdzasumI2(
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
		local_sum = (double*) malloc(dISize(fold));
		created = 1;
	}

	for (i = 0; i < 2 * fold; i++)
		local_sum[i] = 0.0;

	// REGISTER THE MPI TYPE FOR THE EXTENDED INDEXED FP NUMBER
	if (fold != DEFAULT_FOLD)
		dIMPICreate(fold, &myType);
	else
		myType = MPI_IDOUBLE;
   	
	// PEFORM THE LOCAL SUM WITH BLOCK SIZE OF 1024
	dzasumI1(N, x, incx, fold, W, local_sum, local_sum + fold);

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

double prdzasum2(
	MPI_Comm comm, int root,
	int N,
	double complex* x, int incx,
	int fold, int W 
) {
	double *sum = NULL;
	double ret;
	int    me;

	sum = (double*) malloc(dISize(fold));

	prdzasumI2(comm, root, N, x, incx, fold, W, sum, NULL);

	MPI_Comm_rank(comm, &me);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
	ret = Iconv2d1(fold, sum, sum + fold, 1);

	free(sum);
	return ret;
}
