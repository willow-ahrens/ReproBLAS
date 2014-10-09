/*
 *  Created   13/10/25   H.D. Nguyen
 */

#ifndef _REPRODUCIBLE_BLAS_MPI__H_
#define _REPRODUCIBLE_BLAS_MPI__H_

#include <mpi.h>

#include "rblas1.h"
#include "prblas1.h"

//====

void pdcumsum(int me, int procs, int N, double* v, int fold);
void pdcumsum2(int me, int procs, int N, int NB, double* v, int inc, int fold);

void rdsumma1(MPI_Comm comm_row, MPI_Comm comm_col,
		 const int M, const int N, const int K,
		 const int NP, const int DOTNB,
		 int* partitionA,
		 int* partitionB,
		 double alpha,
		 double *A, const int ldA, const int MA, const int NA,
		 double *B, const int ldB, const int MB, const int NB,
		 double beta,
		 double *C, const int ldC, const int MC, const int NC,
		 int fold,
		 double *work);

void rdsumma2(MPI_Comm comm_row, MPI_Comm comm_col,
		 const int M, const int N, const int K,
		 const int NP, const int DOTNB,
		 int* partitionA,
		 int* partitionB,
		 double alpha,
		 double *A, const int ldA, const int MA, const int NA,
		 double *B, const int ldB, const int MB, const int NB,
		 double beta,
		 double *C, const int ldC, const int MC, const int NC,
		 int fold,
		 double *work);
void rdsumma(MPI_Comm comm_row, MPI_Comm comm_col,
	 const int M, const int N, const int K,
	 const int NP,
	 int* partitionA,
	 int* partitionB,
	 double alpha,
	 double *A, const int ldA, const int MA, const int NA,
	 double *B, const int ldB, const int MB, const int NB,
	 double beta,
	 double *C, const int ldC, const int MC, const int NC
);

#endif
