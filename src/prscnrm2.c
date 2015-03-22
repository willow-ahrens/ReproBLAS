/*
 *  Created   13/10/25   H.D. Nguyen
 */

//#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "rblas1.h"
#include "MPI_Indexed.h"
#include "MPI_Indexed/MPI_sIndexed.h"
#include "prblas.h"

void prscnrm2I(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx,
	float* sum
) {
	float BUFFER[1 + sizeof(I_float)/sizeof(float)];
	I_float* local_sum = (I_float*)(BUFFER + 1);

	RMPI_Init();

	sISetZero(*local_sum);

	BUFFER[0] = scnrm2I(N, x, incx, local_sum);

//	sISet(sum+1, local_sum[0]);
	sum[0] = 0;

	if (root >= 0)
		MPI_Reduce(BUFFER,sum,1,MPI_IFLOAT_SCALE, MPI_RNRM2, root, comm);
	else
		MPI_Allreduce(BUFFER,sum, 1,MPI_IFLOAT_SCALE, MPI_RNRM2, comm);
}

float prscnrm2(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx
) {
	float BUFFER[1 + sizeof(I_float)/sizeof(float)];
	I_float* sum = (I_float*)(BUFFER + 1);

	prscnrm2I(comm, root, N, x, incx, BUFFER);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
	return BUFFER[0] * sqrt(Iconv2f(*sum));
}


void prscnrm2I2(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx,
	int fold, int W, float* sum,
	float* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
) {
	int me, i;

	MPI_Datatype myType;

	int created = 0;
	if (local_sum == NULL) {
		local_sum = (float*) malloc(sISize(fold) + sizeof(float));
		created = 1;
	}
	F_CARRY_T* C = (F_CARRY_T*) (local_sum + 1 + fold);

	local_sum[fold] = 0.0;
	for (i = 0; i < fold; i++) {
		local_sum[i] = 0.0;
		C[i] = 0;
	}

	// REGISTER THE MPI TYPE FOR THE EXTENDED INDEXED FP NUMBER
	// REGISTER THE MPI TYPE FOR THE EXTENDED INDEXED FP NUMBER
	if (fold != DEFAULT_FOLD)
		sIMPIScaleCreate(fold, &myType);
	else
		myType = MPI_IFLOAT_SCALE;
   	
	// PEFORM THE LOCAL SUM WITH BLOCK SIZE OF 1024
	local_sum[0] = scnrm2I1(N, x, incx, fold, W, local_sum + 1, C, NULL);

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

float prscnrm22(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx,
	int fold, int W
) {
	float *sum = NULL;
	float *local_sum = NULL;
	float ret;
	int    me;

	sum = (float*) malloc(sizeof(float) + sISize(fold));

	prscnrm2I2(comm, root, N, x, incx, fold, W, sum, NULL);

	MPI_Comm_rank(comm, &me);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
	ret = sum[0] * Iconv2f1(fold, sum + 1, (F_CARRY_T*)(sum + 1 + fold), 1);

	free(sum);
	return ret;
}

