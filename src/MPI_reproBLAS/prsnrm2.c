//#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "MPI_indexed.h"
#include "../MPI_indexed/MPI_sindexed.h"
#include "MPI_reproBLAS.h"
#include "indexedBLAS.h"

void prsnrm2I(
	MPI_Comm comm, int root,
	int N,
	float* x, int incx,
	void* sum
) {
	int i;
	int BUF[3 * DEFAULT_FOLD + 1];
	RMPI_Init();
	void* local_sum = BUF;
	float* psum = (float*) (local_sum);
	float* pcarry = (float*) (psum + 1 + DEFAULT_FOLD);

	for (i = 0; i < DEFAULT_FOLD; i++) {
		psum[i+1] = 0;
		pcarry[i] = 0;
	}

	psum[0] = snrm2I1(N, x, incx, DEFAULT_FOLD, psum + 1, pcarry);

	if (root >= 0)
		MPI_Reduce(local_sum,sum,1,MPI_IFLOAT_SCALE, MPI_RNRM2, root, comm);
	else
		MPI_Allreduce(local_sum,sum, 1,MPI_IFLOAT_SCALE, MPI_RNRM2, comm);
}

float prsnrm2(
	MPI_Comm comm, int root,
	int N,
	float* x, int incx
) {
	// BUFFER
	int sum[16];
	int    me;

	prsnrm2I(comm, root, N, x, incx, sum);

	MPI_Comm_rank(comm, &me);

	float* psum = (float*) (sum);

	if (me == root || root < 0)
		// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
		return psum[0] * sqrt(Iconv2f1(DEFAULT_FOLD,psum + 1, (float*)(psum + 1 + DEFAULT_FOLD), 1));

	return 0.0;
}

void prsnrm2I2(
	MPI_Comm comm, int root,
	int N,
	float* x, int incx,
	int fold, 
	void* sum,
	void* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
) {
	int me, i;

	MPI_Op myOp;
	MPI_Datatype myType;

	int created = 0;
	if (local_sum == NULL) {
		local_sum = malloc(fold * (sizeof(float) + sizeof(float)));
		created = 1;
	}

	float* psum = (float*) (local_sum);
	float* pcarry = (float*) (local_sum + fold+1);

	for (i = 0; i < fold; i++){
		psum[i] = 0.0;
		pcarry[i] = 0;
	}
	psum[fold] = 0.0;

	// REGISTER THE MPI TYPE FOR THE EXTENDED INDEXED FP NUMBER
	if (fold != DEFAULT_FOLD)
		sIMPIScaleCreate(fold, &myType);
	else
		myType = MPI_IFLOAT_SCALE;
   	
	// PEFORM THE LOCAL SUM WITH BLOCK SIZE OF 1024
	psum[0] = snrm2I1(N, x, incx, fold, psum + 1, (float*)(psum + 1 + fold));

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

float prsnrm22(
	MPI_Comm comm, int root,
	int N,
	float* x, int incx,
	int fold
) {
	float ret;
	int    me;

	void* sum = malloc((1+fold) * sizeof(float) + fold * sizeof(float));

	prsnrm2I2(comm, root, N, x, incx, fold, sum, NULL);

	MPI_Comm_rank(comm, &me);

	float* psum = (float*) (sum);

	if (me == root || root < 0)
		// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
		ret = psum[0] * sqrt(Iconv2f1(fold, psum+1, (float*)(psum +1 + fold), 1));
	else
		ret = 0.0;

	free(sum);
	return ret;
}


