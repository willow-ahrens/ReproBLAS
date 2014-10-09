/*
 *  Created   13/10/25   H.D. Nguyen
 */

//#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include "rblas1.h"
#include "IndexedFP/MPIndexedFP.h"
#include "IndexedFP/sIndexedMPI.h"
#include "prblas.h"
#include "types.h"

void prcdotI(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx,
	float complex* y, int incy,
	I_float_Complex* sum,
	int conj
) {
	I_float_Complex local_sum;
	int i;

	RMPI_Init();
	// PEFORM THE LOCAL SUM WITH BLOCK SIZE OF 1024

	if (conj == 1)
		local_sum = cdotcI(N, x, incx, y, incy);
	else
		local_sum = cdotuI(N, x, incx, y, incy);

	cISet(*sum, local_sum);

	// REDUCE THE RESULT TO THE ROOT PROCESSOR
	if (root >= 0)
		MPI_Reduce(&local_sum, sum, 1, MPI_ICOMPLEX,
			MPI_RSUM, root, comm);
	else
		MPI_Allreduce(&local_sum, sum, 1, MPI_ICOMPLEX,
			MPI_RSUM, comm);
}

void prcdotcI(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx,
	float complex* y, int incy,
	I_float_Complex* sum
) {
	prcdotI(comm, root, N, x, incx, y, incy, sum, 1);
}
void prcdotuI(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx,
	float complex* y, int incy,
	I_float_Complex* sum
) {
	prcdotI(comm, root, N, x, incx, y, incy, sum, 0);
}

float complex prcdotc(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx,
	float complex* y, int incy
) {
	I_float_Complex sum;
	int i;

	prcdotI(comm, root, N, x, incx, y, incy, &sum, 1);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
	return Iconv2c(sum);
}

float complex prcdotu(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx,
	float complex* y, int incy
) {
	I_float_Complex sum;
	int i;

	prcdotI(comm, root, N, x, incx, y, incy, &sum, 0);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
	return Iconv2c(sum);
}


void prcdotI2(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx,
	float complex* y, int incy,
	int fold, int W, float complex* sum,
	float complex* local_sum,	// WORKING BUFFER TO STORE LOCAL SUM
	int conj
) {
	int me, i;
	F_CARRY_T* C = (F_CARRY_T*) (local_sum + fold);

	MPI_Datatype myType;
	int created = 0;
	if (local_sum == NULL) {
		local_sum = (float complex*) malloc(cISize(fold));
		created = 1;
	}

	for (i = 0; i < fold; i++) {
		CSET_(local_sum[i], 0.0, 0.0);
		C[2*i] = C[2*i+1] = 0;
	}

	// REGISTER THE MPI TYPE FOR THE EXTENDED INDEXED FP NUMBER
	if (fold != DEFAULT_FOLD)
		sIMPICreate(2*fold, &myType);
	else
		myType = MPI_IDOUBLE;
   	
	// PEFORM THE LOCAL SUM WITH BLOCK SIZE OF 1024
	if (conj == 1)
		cdotcI1(N, x, incx, y, incy, fold, W, local_sum, C);
	else
		cdotuI1(N, x, incx, y, incy, fold, W, local_sum, C);

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

void prcdotcI2(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx,
	float complex* y, int incy,
	int fold, int W, float complex* sum,
	float complex* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
) {
	prcdotI2(comm, root, N, x, incx, y, incy, fold, W, sum, local_sum, 1);
}
void prcdotuI2(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx,
	float complex* y, int incy,
	int fold, int W, float complex* sum,
	float complex* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
) {
	prcdotI2(comm, root, N, x, incx, y, incy, fold, W, sum, local_sum, 0);
}

float complex prcdotc2(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx,
	float complex* y, int incy,
	int fold, int W 
) {
	float complex *sum = NULL;
	float complex ret;
	int    me;

	sum = (float complex*) malloc(cISize(fold));

	prcdotcI2(comm, root, N, x, incx, y, incy, fold, W, sum, NULL);

	MPI_Comm_rank(comm, &me);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
	ret = Iconv2c1(fold, sum, (F_CARRY_T*)(sum + fold), 1);

	free(sum);
	return ret;
}

float complex prcdotu2(
	MPI_Comm comm, int root,
	int N,
	float complex* x, int incx,
	float complex* y, int incy,
	int fold, int W 
) {
	float complex *sum = NULL;
	float complex ret;
	int    me;

	sum = (float complex*) malloc(cISize(fold));

	prcdotuI2(comm, root, N, x, incx, y, incy, fold, W, sum, NULL);

	MPI_Comm_rank(comm, &me);

	// CONVERT THE EXTENDED INDEXED FP NUMBER TO DOUBLE-PRECISION
	ret = Iconv2c1(fold, sum, (F_CARRY_T*)(sum + fold), 1);

	free(sum);
	return ret;
}

