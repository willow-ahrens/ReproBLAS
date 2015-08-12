/*
 *  Created   13/10/25   H.D. Nguyen
 */

//#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "MPI_idxd.h"
#include "../MPI_indexed/MPI_dIndexed.h"
#include "MPI_reproBLAS.h"
#include "indexedBLAS.h"
#include "reproBLAS.h"


void prdgemv(int rank, int nprocs, rblas_order_t Order, rblas_transpose_t TransA, int M, int N, double *myA, int lda, double *myX, int incX, double *Y, int incY){
    (void)lda;
    (void)incX;
    (void)incY;
    double_indexed *myYI;
    double_indexed *YI;
    int i;

    RMPI_Init();
    myYI = (double_indexed*)malloc(M * sizeof(Idouble));
    memset(myYI, 0, M * disize(DEFAULT_FOLD));
    dgemvI(DEFAULT_FOLD, Order, TransA, M, N/nprocs, myA, N/nprocs, myX, 1, myYI, 1);
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
