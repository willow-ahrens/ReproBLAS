/*
 *  Created   13/10/25   H.D. Nguyen
 */

//#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "MPI_indexed.h"
#include "../MPI_indexed/MPI_dIndexed.h"
#include "MPI_reproBLAS.h"
#include "indexedBLAS.h"
#include "reproBLAS.h"
#include "../../tests/common/test_BLAS.h"
#include "../common/common.h"

#define BLK 32

void prbdgemv(int rank, int nprocs, rblas_order_t Order, rblas_transpose_t TransA, int M, int N, double *myA, int lda, double *myX, int incX, double *Y, int incY){
    (void)lda;
    (void)incX;
    (void)incY;
    double_indexed *myYI;
    double_indexed *YI;
    double_indexed *myY;
    int i;
    int ii;

    RMPI_Init();
    myYI = (double_indexed*)malloc(M * sizeof(Idouble));
    myY = (double_indexed*)malloc(M * sizeof(double));
    memset(myY, 0, M * sizeof(double));
    for(i = 0; i < M; i+= BLK){
      CALL_DGEMV(101, 111, MIN((M - i), BLK),N/nprocs,1.0,myA + i * N/nprocs,N/nprocs,myX,incX,1.0,myY + i,incY);

      for(ii = i; ii < M && ii < i + BLK; ii++){
        didconv(DEFAULT_FOLD, myY[ii], myYI + ii * dinum(DEFAULT_FOLD));
      }
    }
    free(myY);
    if(rank == 0){
      YI = (double_indexed*)malloc(M * sizeof(Idouble));
      memset(YI, 0, M * disize(DEFAULT_FOLD));
    }else{
      YI = NULL;
    }

    MPI_Reduce(myYI, YI, M, MPI_IDOUBLE, MPI_RSUM, 0, MPI_COMM_WORLD);
    if(rank == 0){
      for(i = 0; i < M; i++){
        Y[i] = ddiconv(DEFAULT_FOLD, YI + i * dinum(DEFAULT_FOLD));
      }
      free(YI);
    }
    free(myYI);
}
