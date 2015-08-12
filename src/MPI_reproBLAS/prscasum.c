/*
 *  Created   13/10/25   H.D. Nguyen
 */

//#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include "MPI_idxd.h"
#include "../MPI_indexed/MPI_sidxd.h"
#include "MPI_reproBLAS.h"
#include "idxdBLAS.h"

void prscasumI(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx,
	I_float* sum
) {
	I_float local_sum;
	int i;

	RMPI_Init();
	// PEFORM THE LOCAL SUM WITH BLOCK SIZE OF 1024
	local_sum = scasumI(N, x, incx);

	sISet(*sum, local_sum);

	// REDUCE THE RESULT TO THE ROOT PROCESSOR
	if (root >= 0)
		MPI_Reduce(&local_sum, sum, 1, MPI_IFLOAT, MPI_RSUM, root, comm);
	else
		MPI_Allreduce(&local_sum, sum, 1, MPI_IFLOAT, MPI_RSUM, comm);
}

float prscasum(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx
) {
	Ifloat sum;
	int    me;

	prscasumI(comm, root, N, x, incx, &sum);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO FLOAT-PRECISION
	return Iconv2f(sum);
}

void prscasumI2(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx,
	int fold, float* sum,
	float* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
) {
	int me, i;

	MPI_Datatype myType;
	int created = 0;
	if (local_sum == NULL) {
		local_sum = (float*) malloc(sISize(fold));
		created = 1;
	}

	float* C = (float*)(local_sum + fold);
	for (i = 0; i < fold; i++) {
		local_sum[i] = 0;
		C[i] = 0;
	}

	// REGISTER THE MPI TYPE FOR THE EXTENDED INDEXED FP NUMBER
	if (fold != DEFAULT_FOLD)
		sIMPICreate(fold, &myType);
	else
		myType = MPI_IFLOAT;
   	
	// PEFORM THE LOCAL SUM WITH BLOCK SIZE OF 1024
	scasumI1(N, x, incx, fold, local_sum, C, NULL);

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

float prscasum2(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx,
	int fold
) {
	float *sum = NULL;
	float ret;
	int    me;

	sum = (float*) malloc(sISize(fold));

	prscasumI2(comm, root, N, x, incx, fold, sum, NULL);

	MPI_Comm_rank(comm, &me);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO FLOAT-PRECISION
	ret = Iconv2f1(fold, sum, (float*)(sum + fold), 1);

	free(sum);
	return ret;
}
