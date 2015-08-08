#include <indexed.h>
#include <indexedBLAS.h>
#include <reproBLAS.h>

#include "../../config.h"

#include "wrap_daugsum.h"

void wrap_rdgemv(int fold, char Order, char TransA, int M, int N, double alpha, double *A, int lda, double *X, int incX, double beta, double *Y, int incY){
  double_indexed *YI;
  int NY;
  int i;
  if(fold == DIDEFAULTFOLD){
    rdgemv(Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
  }else{
    switch(TransA){
      case 'n':
      case 'N':
        NY = M;
      break;
      default:
        NY = N;
      break;
    }
    YI = (double_indexed*)malloc(NY * incY * disize(fold));
    for(i = 0; i < NY; i++){
      didconv(fold, Y[i * incY] * beta, YI + i * incY * dinum(fold));
    }
    didgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, YI, incY);
    for(i = 0; i < NY; i++){
      Y[i * incY] = ddiconv(fold, YI + i * incY * dinum(fold));
    }
    free(YI);
  }
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
          switch(FillA){
            case util_Mat_Row_Constant:
            case util_Mat_Row_Pos_Big:
            case util_Mat_Row_Pos_Pos_Big:
            case util_Mat_Row_Pos_Neg_Big:
            case util_Mat_Row_Sine:
            case util_Mat_Row_Small_Plus_Increasing_Big:
            case util_Mat_Row_Small_Plus_Rand_Big:
            case util_Mat_Row_N_Cond:
            case util_Mat_Row_Constant_Drop:
            case util_Mat_Row_Sine_Drop:
              switch(Order){
                case 'r':
                case 'R':
                  res[i] += wrap_daugsum_result(M, wrap_daugsum_RDSUM, FillX, A[i] * (RealScaleX*RealAlpha), ImagScaleX, 0, 0.0, 0.0);
                  break;
                default:
                  res[i] += wrap_daugsum_result(M, wrap_daugsum_RDSUM, FillX, A[i * lda] * (RealScaleX*RealAlpha), ImagScaleX, 0, 0.0, 0.0);
                  break;
              }
              break;
            case util_Mat_Row_Pos_Inf:
            case util_Mat_Row_Pos_Pos_Inf:
            case util_Mat_Row_Pos_Neg_Inf:
            case util_Mat_Row_NaN:
            case util_Mat_Row_Pos_Inf_NaN:
            case util_Mat_Row_Pos_Pos_Inf_NaN:
            case util_Mat_Row_Pos_Neg_Inf_NaN:
              for(j = 0; j < M; j++){
                switch(Order){
                  case 'r':
                  case 'R':
                    if(j == i){
                      res[i] += (RealAlpha * X[j * incX]) * A[i];
                    }else{
                      res[i] += (RealAlpha * X[j * incX]) * A[i];
                    }
                    break;
                  default:
                    if(j == i){
                      res[i] += (RealAlpha * X[j * incX]) * A[i * lda];
                    }else{
                      res[i] += (RealAlpha * X[j * incX]) * A[i * lda];
                    }
                    break;
                }
              }
              break;
            default:
              fprintf(stderr, "ReproBLAS error: unknown result for rdgemv(Transpose, A=%s, X=%s)\n", util_mat_fill_descs[FillA], util_vec_fill_descs[FillX]);
              exit(125);
              break;
          }
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
          bound = dibound(fold, N + 1, MAX(fabs(Y[i * incY] * beta), damaxm(N, A + i * lda, 1, alphaX, 1)), res[i * incY]);
          break;
        default:
          bound = dibound(fold, N + 1, MAX(fabs(Y[i * incY] * beta), damaxm(N, A + i, lda, alphaX, 1)), res[i * incY]);
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
          bound = dibound(fold, M + 1, MAX(fabs(Y[i * incY] * beta), damaxm(M, A + i, lda, alphaX, 1)), res[i * incY]);
          break;
        default:
          bound = dibound(fold, M + 1, MAX(fabs(Y[i * incY] * beta), damaxm(M, A + i * lda, 1, alphaX, 1)), res[i * incY]);
          break;
      }
      break;
  }
  free(alphaX);
  return bound;
}
