#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <indexed.h>
#include <indexedBLAS.h>
#include <reproBLAS.h>

#include "../common/test_opt.h"
#include "../common/test_matvec_fill_header.h"

#include "wrap_rdgemv.h"

static opt_option max_blocks;
static opt_option shuffles;
static opt_option fold;

static void verify_rdgemv_options_initialize(void){
  max_blocks._int.header.type       = opt_int;
  max_blocks._int.header.short_name = 'B';
  max_blocks._int.header.long_name  = "blocks";
  max_blocks._int.header.help       = "maximum number of blocks";
  max_blocks._int.required          = 0;
  max_blocks._int.min               = 1;
  max_blocks._int.max               = INT_MAX;
  max_blocks._int.value             = 1024;

  shuffles._int.header.type       = opt_int;
  shuffles._int.header.short_name = 'S';
  shuffles._int.header.long_name  = "shuffles";
  shuffles._int.header.help       = "number of times to shuffle";
  shuffles._int.required          = 0;
  shuffles._int.min               = 0;
  shuffles._int.max               = INT_MAX;
  shuffles._int.value             = 5;

  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 2;
  fold._int.max               = DIMAXFOLD;
  fold._int.value             = DIDEFAULTFOLD;
}

int verify_dgemv_reproducibility(int fold, char Order, char TransA, int M, int N, int NX, int NY, double alpha, double *A, int lda, double* X, int incX, double beta, double *Y, double_indexed *YI, int incY, double *ref, double_indexed *Iref, int max_num_blocks) {

  int i;
  int num_blocks = 1;
  int block_N;

  double *res = malloc(NY * incY * sizeof(double));
  double_indexed *Ires = malloc(NY * incY * disize(fold));

  num_blocks = 1;
  while (num_blocks < N && num_blocks <= max_num_blocks) {
    memcpy(res, Y, NY * incY * sizeof(double));
    memcpy(Ires, YI, NY * incY * disize(fold));
    if (num_blocks == 1){
      wrap_rdgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, beta, res, incY);
    }else {
      switch(TransA){
        case 'n':
        case 'N':
          {
            block_N = (N + num_blocks - 1) / num_blocks;
            for (i = 0; i < N; i += block_N) {
              block_N = block_N < N - i ? block_N : (N-i);
              switch(Order){
                case 'r':
                case 'R':
                  didgemv(fold, Order, TransA, M, block_N, alpha, A + i, lda, X + i * incX, incX, Ires, incY);
                  break;
                default:
                  didgemv(fold, Order, TransA, M, block_N, alpha, A + i * lda, lda, X + i * incX, incX, Ires, incY);
                  break;
              }
            }
          }
          break;
        default:
          {
            block_N = (M + num_blocks - 1) / num_blocks;
            for (i = 0; i < M; i += block_N) {
              block_N = block_N < M - i ? block_N : (M-i);
              switch(Order){
                case 'r':
                case 'R':
                  didgemv(fold, Order, TransA, block_N, N, alpha, A + i * lda, lda, X + i * incX, incX, Ires, incY);
                  break;
                default:
                  didgemv(fold, Order, TransA, block_N, N, alpha, A + i, lda, X + i * incX, incX, Ires, incY);
                  break;
              }
            }
          }
          break;
      }
      for(i = 0; i < NY; i++){
        res[i * incY] = ddiconv(fold, Ires + i * incY * dinum(fold));
      }
    }
    for(i = 0; i < NY; i++){
      if(res[i * incY] != ref[i * incY]){
        printf("rdgemv(A, X, Y)[num_blocks=%d,block_N=%d] = %g != %g\n", num_blocks, block_N, res[i * incY], ref[i * incY]);
        if (num_blocks != 1) {
          printf("Ref I_double:\n");
          diprint(fold, Iref + i * incY * dinum(fold));
          printf("\nRes I_double:\n");
          diprint(fold, Ires + i * incY * dinum(fold));
          printf("\n");
        }
        return 1;
      }
    }
    num_blocks *= 2;
  }
  return 0;
}

int matvec_fill_show_help(void){
  verify_rdgemv_options_initialize();

  opt_show_option(fold);
  opt_show_option(max_blocks);
  opt_show_option(shuffles);
  return 0;
}

const char* matvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  verify_rdgemv_options_initialize();

  opt_eval_option(argc, argv, &fold);
  opt_eval_option(argc, argv, &max_blocks);

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Verify rdgemv reproducibility fold=%d", fold._int.value);
  return name_buffer;
}

int matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillX, double RealScaleX, double ImagScaleX, int incX, double RealBeta, double ImagBeta, int FillY, double RealScaleY, double ImagScaleY, int incY){
  (void)ImagAlpha;
  (void)ImagBeta;
  int rc = 0;
  int i;

  verify_rdgemv_options_initialize();

  opt_eval_option(argc, argv, &fold);
  opt_eval_option(argc, argv, &max_blocks);
  opt_eval_option(argc, argv, &shuffles);

  util_random_seed();
  int NX;
  int NY;
  char NTransA;
  switch(TransA){
    case 'n':
    case 'N':
      NX = N;
      NY = M;
      NTransA = 'T';
    break;
    default:
      NX = M;
      NY = N;
      NTransA = 'N';
    break;
  }

  double *A  = util_dmat_alloc(Order, M, N, lda);
  double *X  = util_dvec_alloc(NX, incX);
  double *Y  = util_dvec_alloc(NY, incY);
  double_indexed *YI = (double_indexed*)malloc(NY * incY * disize(fold._int.value));

  int *P;

  util_dmat_fill(Order, 'n', M, N, A, lda, FillA, RealScaleA, ImagScaleA);
  util_dvec_fill(NX, X, incX, FillX, RealScaleX, ImagScaleX);
  util_dvec_fill(NY, Y, incY, FillY, RealScaleY, ImagScaleY);
  for(i = 0; i < NY; i++){
    didconv(fold._int.value, Y[i * incY] * RealBeta, YI + i * incY * dinum(fold._int.value));
  }
  double *ref  = (double*)malloc(NY * incY * sizeof(double));
  double_indexed *Iref = (double_indexed*)malloc(NY * incY * disize(fold._int.value));

  //compute with unpermuted data
  memcpy(ref, Y, NY * incY * sizeof(double));
  memcpy(Iref, YI, NY * incY * disize(fold._int.value));

  wrap_rdgemv(fold._int.value, Order, TransA, M, N, RealAlpha, A, lda, X, incX, RealBeta, ref, incY);
  didgemv(fold._int.value, Order, TransA, M, N, RealAlpha, A, lda, X, incX, Iref, incY);

  P = util_identity_permutation(NX);
  util_dvec_reverse(NX, X, incX, P, 1);
  util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_dgemv_reproducibility(fold._int.value, Order, TransA, M, N, NX, NY, RealAlpha, A, lda, X, incX, RealBeta, Y, YI, incY, ref, Iref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_dvec_sort(NX, X, incX, P, 1, util_Increasing);
  util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_dgemv_reproducibility(fold._int.value, Order, TransA, M, N, NX, NY, RealAlpha, A, lda, X, incX, RealBeta, Y, YI, incY, ref, Iref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_dvec_sort(NX, X, incX, P, 1, util_Decreasing);
  util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_dgemv_reproducibility(fold._int.value, Order, TransA, M, N, NX, NY, RealAlpha, A, lda, X, incX, RealBeta, Y, YI, incY, ref, Iref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_dvec_sort(NX, X, incX, P, 1, util_Increasing_Magnitude);
  util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_dgemv_reproducibility(fold._int.value, Order, TransA, M, N, NX, NY, RealAlpha, A, lda, X, incX, RealBeta, Y, YI, incY, ref, Iref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_dvec_sort(NX, X, incX, P, 1, util_Decreasing_Magnitude);
  util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_dgemv_reproducibility(fold._int.value, Order, TransA, M, N, NX, NY, RealAlpha, A, lda, X, incX, RealBeta, Y, YI, incY, ref, Iref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  for(i = 0; i < shuffles._int.value; i++){
    P = util_identity_permutation(NX);
    util_dvec_shuffle(NX, X, incX, P, 1);
    util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
    free(P);

    rc = verify_dgemv_reproducibility(fold._int.value, Order, TransA, M, N, NX, NY, RealAlpha, A, lda, X, incX, RealBeta, Y, YI, incY, ref, Iref, max_blocks._int.value);
    if(rc != 0){
      return rc;
    }
  }

  free(A);
  free(X);
  free(Y);
  free(ref);
  free(Iref);

  return rc;
}
