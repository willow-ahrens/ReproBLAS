#include <idxd.h>
#include <indexedBLAS.h>
#include <reproBLAS.h>

#include "../../config.h"

#include "wrap_daugsum.h"

void wrap_rdgemv(int fold, char Order, char TransA, int M, int N, double alpha, double *A, int lda, double *X, int incX, double beta, double *Y, int incY){
  double_indexed *YI;
  int opM;
  int i;
  if(fold == DIDEFAULTFOLD){
    rdgemv(Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
  }else{
    switch(TransA){
      case 'n':
      case 'N':
        opM = M;
      break;
      default:
        opM = N;
      break;
    }
    YI = (double_indexed*)malloc(opM * incY * idxd_disize(fold));
    for(i = 0; i < opM; i++){
      idxd_didconv(fold, Y[i * incY] * beta, YI + i * incY * idxd_dinum(fold));
    }
    didgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, YI, incY);
    for(i = 0; i < opM; i++){
      Y[i * incY] = idxd_ddiconv(fold, YI + i * incY * idxd_dinum(fold));
    }
    free(YI);
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
  for(i = 0; i < opN; i++){
    opX[i] = alpha * X[i * incX];
  }
  for(i = 0; i < opM; i++){
    if(beta == 0.0){
      idxd_disetzero(fold, YI);
    }else{
      idxd_didconv(fold, Y[i * incY] * beta, YI);
    }
    if(alpha != 0.0){
      switch(Order){
        case 'r':
        case 'R':
          diddot(fold, opN, opA + i * opN, 1, opX, 1, YI);
          break;
        default:
          diddot(fold, opN, opA + i, opN, opX, 1, YI);
          break;
      }
    }
    Y[i * incY] = idxd_ddiconv(fold, YI);
  }
  free(YI);
  free(opA);
}

double* wrap_rdgemv_result(char Order, char TransA, int M, int N, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, double *A, int lda, int FillX, double RealScaleX, double ImagScaleX, double *X, int incX, double RealBeta, double ImagBeta, double *Y, int incY){
  int i;
  int j;
  double *res;
  switch(TransA){
    case 'n':
    case 'N':
      res = (double*)malloc(M * sizeof(double));
      for(i = 0; i < M; i++){
        res[i] = Y[i * incY] * RealBeta;
        if(FillA == util_Mat_Identity){
          for(j = 0; j < N; j++){
            if(j == i){
              res[i] += (RealAlpha * X[j * incX]) * 1.0;
            }else{
              res[i] += (RealAlpha * X[j * incX]) * 0.0;
            }
          }
        }else{
          res[i] += wrap_daugsum_result(N, wrap_daugsum_RDDOT, FillX, RealScaleX*RealAlpha, ImagScaleX, (util_vec_fill_t)FillA, RealScaleA, ImagScaleA);
        }
      }
      break;
    default:
      res = (double*)malloc(N * sizeof(double));
      for(i = 0; i < N; i++){
        res[i] = Y[i * incY] * RealBeta;
        if(FillA == util_Mat_Identity){
          for(j = 0; j < M; j++){
            if(j == i){
              res[i] += (RealAlpha * X[j * incX]) * 1.0;
            }else{
              res[i] += (RealAlpha * X[j * incX]) * 0.0;
            }
          }
        }else{
          res[i] += wrap_daugsum_result(M, wrap_daugsum_RDDOT, FillX, RealScaleX*RealAlpha, ImagScaleX, (util_vec_fill_t)FillA, RealScaleA, ImagScaleA);
        }
      }
      break;
  }
  return res;
}

double wrap_rdgemv_bound(int fold, char Order, char TransA, int M, int N, double alpha, double *A, int lda, double *X, int incX, double beta, double *Y, int incY, double *res, double *ref, int i){
  int j;
  double *alphaX;
  double bound;
  switch(TransA){
    case 'n':
    case 'N':
      alphaX = (double*)malloc(N * sizeof(double));
      for(j = 0; j < N; j++){
        alphaX[j] = X[j * incX] * alpha;
      }
      switch(Order){
        case 'r':
        case 'R':
          bound = idxd_dibound(fold, N + 1, MAX(fabs(Y[i * incY] * beta), damaxm(N, A + i * lda, 1, alphaX, 1)), res[i * incY]);
          break;
        default:
          bound = idxd_dibound(fold, N + 1, MAX(fabs(Y[i * incY] * beta), damaxm(N, A + i, lda, alphaX, 1)), res[i * incY]);
          break;
      }
      break;
    default:
      alphaX = (double*)malloc(M * sizeof(double));
      for(j = 0; j < M; j++){
        alphaX[j] = X[j * incX] * ((double*)&alpha)[0];
      }
      switch(Order){
        case 'r':
        case 'R':
          bound = idxd_dibound(fold, M + 1, MAX(fabs(Y[i * incY] * beta), damaxm(M, A + i, lda, alphaX, 1)), res[i * incY]);
          break;
        default:
          bound = idxd_dibound(fold, M + 1, MAX(fabs(Y[i * incY] * beta), damaxm(M, A + i * lda, 1, alphaX, 1)), res[i * incY]);
          break;
      }
      break;
  }
  free(alphaX);
  return bound;
}
