//#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include "rblas1.h"
#include "IndexedFP/sIndexedMPI.h"
#include "IndexedFP/MPIndexedFP.h"
#include "prblas.h"

void prsdotI(
	MPI_Comm comm, int root,
	int N,
	float* x, int incx,
	float* y, int incy,
	I_float* sum
) {

	RMPI_Init();
	Ifloat local_sum;

	// PEFORM THE LOCAL COMPUTATION
	local_sum = sdotI(N, x, incx, y, incy);
	sISet(*sum, local_sum);

	// REDUCE THE RESULT TO THE ROOT PROCESSOR
	if (root >= 0)
		MPI_Reduce(&local_sum, sum, 1, MPI_IFLOAT, MPI_RSUM, root, comm);
	else
		MPI_Allreduce(&local_sum, sum, 1, MPI_IFLOAT, MPI_RSUM, comm);
}

float prsdot(
	MPI_Comm comm, int root,
	int N,
	float* x, int incx,
	float* y, int incy
) {
	// BUFFER
	Ifloat sum;
	int    me;

	prsdotI(comm, root, N, x, incx, y, incy, &sum);

	MPI_Comm_rank(comm, &me);

	if (me == root || root < 0)
		// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
		return Iconv2f(sum);

	return 0.0;
}

void prsdotI2(
	MPI_Comm comm, int root,
	int N,
	float* x, int incx,
	float* y, int incy,
	int fold, int W,
	void* sum,
	void* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
) {
	int me, i;

	MPI_Op myOp;
	MPI_Datatype myType;

	int created = 0;
	if (local_sum == NULL) {
		local_sum = malloc(fold * (sizeof(float) + sizeof(F_CARRY_T)));
		created = 1;
	}

	float* psum = (float*) (local_sum);
	F_CARRY_T* pcarry = (F_CARRY_T*) (psum + fold);

	for (i = 0; i < fold; i++){
		psum[i] = 0.0;
		pcarry[i] = 0;
	}

	// REGISTER THE MPI TYPE FOR THE EXTENDED INDEXED FP NUMBER
	if (fold != DEFAULT_FOLD)
		sIMPICreate(fold, &myType);
	else
		myType = MPI_IFLOAT;
   	
	// PEFORM THE LOCAL SUM WITH BLOCK SIZE OF 1024
	sdotI1(N, x, incx, y, incy, fold, W, psum, (F_CARRY_T*)(psum + fold));

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

float prsdot2(
	MPI_Comm comm, int root,
	int N,
	float* x, int incx,
	float* y, int incy,
	int fold, int W 
) {
	float ret;
	int    me;

	void* sum = malloc(fold * (sizeof(float) + sizeof(F_CARRY_T)));

	prsdotI2(comm, root, N, x, incx, y, incy, fold, W, sum, NULL);

	MPI_Comm_rank(comm, &me);

	float* psum = (float*) (sum);

	if (me == root || root < 0)
		// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
		ret = Iconv2f1(fold, psum, (F_CARRY_T*)(psum + fold), 1);
	else
		ret = 0.0;

	free(sum);
	return ret;
}

