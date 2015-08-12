#include <idxd.h>
#include <indexedBLAS.h>
#include <reproBLAS.h>

#include "../../config.h"

#include "wrap_zaugsum.h"

void wrap_rzgemv(int fold, char Order, char TransA, int M, int N, double complex *alpha, double complex *A, int lda, double complex *X, int incX, double complex *beta, double complex *Y, int incY){
  double_indexed *YI;
  double complex betaY;
  int NY;
  int i;
  if(fold == DIDEFAULTFOLD){
    rzgemv(Order, TransA, M, N, (void*)alpha, (void*)A, lda, (void*)X, incX, (void*)beta, (void*)Y, incY);
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
    YI = (double_indexed*)malloc(NY * incY * zisize(fold));
    for(i = 0; i < NY; i++){
      betaY = Y[i * incY] * (*beta);
      zizconv(fold, &betaY, YI + i * incY * zinum(fold));
    }
    zizgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, YI, incY);
    for(i = 0; i < NY; i++){
      zziconv_sub(fold, YI + i * incY * zinum(fold), Y + i * incY);
    }
    free(YI);
  }
}

double complex* wrap_rzgemv_result(char Order, char TransA, int M, int N, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, double complex *A, int lda, int FillX, double RealScaleX, double ImagScaleX, double complex *X, int incX, double RealBeta, double ImagBeta, double complex *Y, int incY){
  int i;
  int j;
  double complex *res;
  switch(TransA){
    case 'n':
    case 'N':
      res = (double complex*)malloc(M * sizeof(double complex));
      for(i = 0; i < M; i++){
        res[i] = Y[i * incY] * (RealBeta + I * ImagBeta);
        if(FillA == util_Mat_Identity){
          for(j = 0; j < N; j++){
            if(j == i){
              res[i] += ((RealAlpha + I * ImagAlpha)* X[j * incX]) * 1.0;
            }else{
              res[i] += ((RealAlpha + I * ImagAlpha) * X[j * incX]) * 0.0;
            }
          }
        }else{
          res[i] += wrap_zaugsum_result(N, wrap_zaugsum_RZDOTU, FillX, RealScaleX * RealAlpha - ImagScaleX * ImagAlpha, ImagScaleX * RealAlpha + RealScaleX * ImagAlpha, (util_vec_fill_t)FillA, RealScaleA, ImagScaleA);
        }
      }
      break;
    default:
      res = (double complex*)malloc(N * sizeof(double complex));
      for(i = 0; i < N; i++){
        res[i] = Y[i * incY] * (RealBeta + I * ImagBeta);
        if(FillA == util_Mat_Identity){
          for(j = 0; j < M; j++){
            if(j == i){
              res[i] += ((RealAlpha + I * ImagAlpha) * X[j * incX]) * 1.0;
            }else{
              res[i] += ((RealAlpha + I * ImagAlpha) * X[j * incX]) * 0.0;
            }
          }
        }else{
          res[i] += wrap_zaugsum_result(M, wrap_zaugsum_RZDOTU, FillX, RealScaleX * RealAlpha - ImagScaleX * ImagAlpha, ImagScaleX * RealAlpha + RealScaleX * ImagAlpha, (util_vec_fill_t)FillA, RealScaleA, ImagScaleA);
        }
      }
      break;
  }
  return res;
}

double complex wrap_rzgemv_bound(int fold, char Order, char TransA, int M, int N, double complex *alpha, double complex *A, int lda, double complex *X, int incX, double complex *beta, double complex *Y, int incY, double complex *res, double complex *ref, int i){
  int j;
  double complex *alphaX;
  double complex betaY;
  double complex bound;
  double complex amaxm;
  double *bound_base = (double*)&bound;
  switch(TransA){
    case 'n':
    case 'N':
      alphaX = (double complex*)malloc(N * sizeof(double complex));
      for(j = 0; j < N; j++){
        alphaX[j] = X[j * incX] * (*alpha);
      }
      switch(Order){
        case 'r':
        case 'R':
          zamaxm_sub(N, A + i * lda, 1, alphaX, 1, &amaxm);
          break;
        default:
          zamaxm_sub(N, A + i, lda, alphaX, 1, &amaxm);
          break;
      }
      betaY = Y[i * incY] * (*beta);
      bound_base[0] = dibound(fold, N + 1, MAX(fabs(creal(betaY)), creal(amaxm)), creal(res[i * incY]));
      bound_base[1] = dibound(fold, N + 1, MAX(fabs(cimag(betaY)), cimag(amaxm)), cimag(res[i * incY]));
      break;
    default:
      alphaX = (double complex*)malloc(M * sizeof(double complex));
      for(j = 0; j < M; j++){
        alphaX[j] = X[j * incX] * (*alpha);
      }
      switch(Order){
        case 'r':
        case 'R':
          zamaxm_sub(M, A + i, lda, alphaX, 1, &amaxm);
          break;
        default:
          zamaxm_sub(M, A + i * lda, 1, alphaX, 1, &amaxm);
          break;
      }
      betaY = Y[i * incY] * (*beta);
      bound_base[0] = dibound(fold, M + 1, MAX(fabs(creal(betaY)), creal(amaxm)), creal(res[i * incY]));
      bound_base[1] = dibound(fold, M + 1, MAX(fabs(cimag(betaY)), cimag(amaxm)), cimag(res[i * incY]));
      break;
  }
  free(alphaX);
  return bound;
}
