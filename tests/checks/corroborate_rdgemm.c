#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <idxd.h>
#include <indexedBLAS.h>
#include <reproBLAS.h>

#include "../common/test_opt.h"
#include "../common/test_matmat_fill_header.h"

#include "wrap_rdgemm.h"

static opt_option max_blocks;
static opt_option shuffles;
static opt_option fold;

static void corroborate_rdgemm_options_initialize(void){
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

int corroborate_rdgemm(int fold, char Order, char TransA, char TransB, int M, int N, int K, double alpha, double *A, int lda, double* B, int ldb, double beta, double *C, double_indexed *CI, int ldc, double *ref, int max_num_blocks) {

  int i;
  int j;
  int k;
  int num_blocks = 1;
  int block_K;

  double *res;
  double_indexed *Ires;
  double *tmpA;
  double *tmpB;
  int CNM;

  switch(Order){
    case 'r':
    case 'R':
      CNM = M * ldc;
      break;
    default:
      CNM = ldc * N;
      break;
  }
  res = malloc(CNM * sizeof(double));
  Ires = malloc(CNM * idxd_disize(fold));

  num_blocks = 1;
  while (num_blocks < K && num_blocks <= max_num_blocks) {
    memcpy(res, C, CNM * sizeof(double));
    memcpy(Ires, CI, CNM * idxd_disize(fold));
    if (num_blocks == 1){
      wrap_rdgemm(fold, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, res, ldc);
    }else {
      block_K = (K + num_blocks - 1) / num_blocks;
      for (k = 0; k < K; k += block_K) {
        block_K = block_K < K - k ? block_K : (K-k);
        switch(Order){
          case 'r':
          case 'R':
            switch(TransA){
              case 'n':
              case 'N':
                tmpA = A + k;
                break;
              default:
                tmpA = A + k * lda;
                break;
            }
            switch(TransB){
              case 'n':
              case 'N':
                tmpB = B + k * ldb;
                break;
              default:
                tmpB = B + k;
                break;
            }
            break;
          default:
            switch(TransA){
              case 'n':
              case 'N':
                tmpA = A + k * lda;
                break;
              default:
                tmpA = A + k;
                break;
            }
            switch(TransB){
              case 'n':
              case 'N':
                tmpB = B + k;
                break;
              default:
                tmpB = B + k * ldb;
                break;
            }
            break;
        }
        didgemm(fold, Order, TransA, TransB, M, N, block_K, alpha, tmpA, lda, tmpB, ldb, Ires, ldc);
      }
      for(i = 0; i < M; i++){
        for(j = 0; j < N; j++){
          switch(Order){
            case 'r':
            case 'R':
              res[i * ldc + j] = idxd_ddiconv(fold, Ires + (i * ldc + j) * idxd_dinum(fold));
              break;
            default:
              res[i + j * ldc] = idxd_ddiconv(fold, Ires + (i + j * ldc) * idxd_dinum(fold));
              break;
          }
        }
      }
    }
    for(i = 0; i < M; i++){
      for(j = 0; j < N; j++){
        switch(Order){
          case 'r':
          case 'R':
            if(res[i * ldc + j] != ref[i * ldc + j]){
              printf("rdgemm(A, X, Y)[num_blocks=%d,block_K=%d] = %g != %g\n", num_blocks, block_K, res[i * ldc + j], ref[i * ldc + j]);
              return 1;
            }
            break;
          default:
            if(res[i + j * ldc] != ref[i + j * ldc]){
              printf("rdgemm(A, X, Y)[num_blocks=%d,block_K=%d] = %g != %g\n", num_blocks, block_K, res[i + j * ldc], ref[i + j * ldc]);
              return 1;
            }
            break;
        }
      }
    }
    num_blocks *= 2;
  }
  return 0;
}

int matmat_fill_show_help(void){
  corroborate_rdgemm_options_initialize();

  opt_show_option(fold);
  opt_show_option(max_blocks);
  opt_show_option(shuffles);
  return 0;
}

const char* matmat_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  corroborate_rdgemm_options_initialize();

  opt_eval_option(argc, argv, &fold);
  opt_eval_option(argc, argv, &max_blocks);

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Corroborate rdgemm fold=%d", fold._int.value);
  return name_buffer;
}

int matmat_fill_test(int argc, char** argv, char Order, char TransA, char TransB, int M, int N, int K, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillB, double RealScaleB, double ImagScaleB, int ldb, double RealBeta, double ImagBeta, int FillC, double RealScaleC, double ImagScaleC, int ldc){
  (void)ImagAlpha;
  (void)ImagBeta;
  int rc = 0;
  int i;
  int j;

  corroborate_rdgemm_options_initialize();

  opt_eval_option(argc, argv, &fold);
  opt_eval_option(argc, argv, &max_blocks);
  opt_eval_option(argc, argv, &shuffles);

  util_random_seed();
  char NTransA;
  switch(TransA){
    case 'n':
    case 'N':
      NTransA = 'T';
    break;
    default:
      NTransA = 'N';
    break;
  }

  double *A  = util_dmat_alloc(Order, M, K, lda);
  double *B  = util_dmat_alloc(Order, K, N, ldb);
  double *C  = util_dmat_alloc(Order, M, N, ldc);
  int CNM;
  switch(Order){
    case 'r':
    case 'R':
      CNM = M * ldc;
      break;
    default:
      CNM = ldc * N;
      break;
  }
  double_indexed *CI = malloc(CNM * idxd_disize(fold._int.value));

  int *P;

  util_dmat_fill(Order, NTransA, M, K, A, lda, FillA, RealScaleA, ImagScaleA);
  util_dmat_fill(Order, TransB, K, N, B, ldb, FillB, RealScaleB, ImagScaleB);
  util_dmat_fill(Order, 'n', M, N, C, ldc, FillC, RealScaleC, ImagScaleC);
  for(i = 0; i < M; i++){
    for(j = 0; j < N; j++){
      switch(Order){
        case 'r':
        case 'R':
          idxd_didconv(fold._int.value, C[i * ldc + j] * RealBeta, CI + (i * ldc + j) * idxd_dinum(fold._int.value));
          break;
        default:
          idxd_didconv(fold._int.value, C[i + j * ldc] * RealBeta, CI + (i + j * ldc) * idxd_dinum(fold._int.value));
          break;
      }
    }
  }
  double *ref  = (double*)malloc(CNM * sizeof(double));

  //compute with unpermuted data
  memcpy(ref, C, CNM * sizeof(double));

  wrap_ref_rdgemm(fold._int.value, Order, TransA, TransB, M, N, K, RealAlpha, A, lda, B, ldb, RealBeta, ref, ldc);

  rc = corroborate_rdgemm(fold._int.value, Order, TransA, TransB, M, N, K, RealAlpha, A, lda, B, ldb, RealBeta, C, CI, ldc, ref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(K);
  util_dmat_row_reverse(Order, NTransA, M, K, A, lda, P, 1);
  util_dmat_row_permute(Order, TransB, K, N, B, ldb, P, 1, NULL, 1);
  free(P);

  rc = corroborate_rdgemm(fold._int.value, Order, TransA, TransB, M, N, K, RealAlpha, A, lda, B, ldb, RealBeta, C, CI, ldc, ref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(K);
  util_dmat_row_sort(Order, NTransA, M, K, A, lda, P, 1, util_Increasing, 0);
  util_dmat_row_permute(Order, TransB, K, N, B, ldb, P, 1, NULL, 1);
  free(P);

  rc = corroborate_rdgemm(fold._int.value, Order, TransA, TransB, M, N, K, RealAlpha, A, lda, B, ldb, RealBeta, C, CI, ldc, ref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(K);
  util_dmat_row_sort(Order, NTransA, M, K, A, lda, P, 1, util_Decreasing, 0);
  util_dmat_row_permute(Order, TransB, K, N, B, ldb, P, 1, NULL, 1);
  free(P);

  rc = corroborate_rdgemm(fold._int.value, Order, TransA, TransB, M, N, K, RealAlpha, A, lda, B, ldb, RealBeta, C, CI, ldc, ref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(K);
  util_dmat_row_sort(Order, NTransA, M, K, A, lda, P, 1, util_Increasing_Magnitude, 0);
  util_dmat_row_permute(Order, TransB, K, N, B, ldb, P, 1, NULL, 1);
  free(P);

  rc = corroborate_rdgemm(fold._int.value, Order, TransA, TransB, M, N, K, RealAlpha, A, lda, B, ldb, RealBeta, C, CI, ldc, ref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(K);
  util_dmat_row_sort(Order, NTransA, M, K, A, lda, P, 1, util_Decreasing_Magnitude, 0);
  util_dmat_row_permute(Order, TransB, K, N, B, ldb, P, 1, NULL, 1);
  free(P);

  rc = corroborate_rdgemm(fold._int.value, Order, TransA, TransB, M, N, K, RealAlpha, A, lda, B, ldb, RealBeta, C, CI, ldc, ref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  for(i = 0; i < shuffles._int.value; i++){
    P = util_identity_permutation(K);
    util_dmat_row_shuffle(Order, NTransA, M, K, A, lda, P, 1);
    util_dmat_row_permute(Order, TransB, K, N, B, ldb, P, 1, NULL, 1);
    free(P);

    rc = corroborate_rdgemm(fold._int.value, Order, TransA, TransB, M, N, K, RealAlpha, A, lda, B, ldb, RealBeta, C, CI, ldc, ref, max_blocks._int.value);
    if(rc != 0){
      return rc;
    }
  }

  free(A);
  free(B);
  free(C);
  free(CI);
  free(ref);

  return rc;
}
