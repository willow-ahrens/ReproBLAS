/*
 *  Created   13/10/25   H.D. Nguyen
 */

//#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include "MPI_indexed.h"
#include "../MPI_indexed/MPI_dindexed.h"
#include "MPI_reproBLAS.h"
#include "indexedBLAS.h"
#include "reproBLAS.h"


void prdgemv(int px, int npx, int py, int npy, rblas_order_t Order, rblas_transpose_t TransA, int M, int N, double *myA, int lda, double *myX, int incX, double *myY, int incY){
    Idouble *myYI;
    Idouble *YI;
    int i;

    RMPI_Init();
    myYI = (Idouble*)malloc(M * sizeof(Idouble));
    for(i = 0; i < M; i++){
      sISetZero(myYI[i]);
    }
    dgemvI(Order, TransA, M/npy, N/npx, myA, N/npx, myX, 1, YI, 1, DEFAULT_FOLD);

    if(py == 0){
      YI = (Idouble*)malloc(M * sizeof(Idouble));
      for(i = 0; i < M; i++){
        sISetZero(YI[i]);
      }
    }

    MPI_Reduce(myYI, YI, M, MPI_IDOUBLE, MPI_RSUM, 0, MPI_COMM_WORLD);

    if(py == 0){
      for(i = 0; i < M; i++){
        myY[i] = ddiconv(YI + i, DEFAULT_FOLD);
      }
    }
}
