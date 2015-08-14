#include <idxd.h>
#include <idxdBLAS.h>

#include "../common/test_util.h"

#include "../../config.h"

void wrap_rdgemv(int fold, char Order, char TransA, int M, int N, double alpha, double *A, int lda, double *X, int incX, double beta, double *Y, int incY){
  if(fold == DIDEFAULTFOLD){
    reproBLAS_dgemv(Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
  }else{
    reproBLAS_rdgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
  }
}

void wrap_ref_rdgemv(int fold, char Order, char TransA, int M, int N, double alpha, double *A, int lda, double *X, int incX, double beta, double *Y, int incY){
  int opM;
  int opN;
  double *opA;
  double *opX;
  double_indexed *YI;
  int i;
  switch(TransA){
    case 'n':
    case 'N':
      opM = M;
      opN = N;
      break;
    default:
      opM = N;
      opN = M;
      break;
  }
  opA = util_dmat_op(Order, TransA, opM, opN, A, lda);
  YI = idxd_dialloc(fold);
  opX = (double*)malloc(opN * sizeof(double));
  for(i = 0; i < opM; i++){
    if(beta == 0.0){
      idxd_disetzero(fold, YI);
    }else{
      idxd_didconv(fold, Y[i * incY] * beta, YI);
    }
    if(alpha != 0.0){
      for(i = 0; i < opN; i++){
        opX[i] = alpha * X[i * incX];
      }
      switch(Order){
        case 'r':
        case 'R':
          idxdBLAS_diddot(fold, opN, opA + i * opN, 1, opX, 1, YI);
          break;
        default:
          idxdBLAS_diddot(fold, opN, opA + i, opN, opX, 1, YI);
          break;
      }
    }
    Y[i * incY] = idxd_ddiconv(fold, YI);
  }
  free(YI);
  free(opA);
}
