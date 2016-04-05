#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <idxd.h>
#include <idxdBLAS.h>
#include <reproBLAS.h>

#include "../common/test_opt.h"
#include "../common/test_matmat_fill_header.h"

#include "wrap_rzgemm.h"

static opt_option max_blocks;
static opt_option shuffles;
static opt_option fold;

static void corroborate_rzgemm_options_initialize(void){
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
  fold._int.max               = idxd_DIMAXFOLD;
  fold._int.value             = DIDEFAULTFOLD;
}

int corroborate_rzgemm(int fold, char Order, char TransA, char TransB, int M, int N, int K, double complex *alpha, double complex *A, int lda, double complex *B, int ldb, double complex *beta, double complex *C, double_complex_indexed *CI, int ldc, double complex *ref, int max_num_blocks) {

  int i;
  int j;
  int k;
  int num_blocks = 1;
  int block_K;

  double complex *res;
  double_complex_indexed *Ires;
  double complex *tmpA;
  double complex *tmpB;
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
  res = malloc(CNM * sizeof(double complex));
  Ires = malloc(CNM * idxd_zisize(fold));

  num_blocks = 1;
  while (num_blocks < K && num_blocks <= max_num_blocks) {
    memcpy(res, C, CNM * sizeof(double complex));
    memcpy(Ires, CI, CNM * idxd_zisize(fold));
    if (num_blocks == 1){
      wrap_rzgemm(fold, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, res, ldc);
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
        idxdBLAS_zizgemm(fold, Order, TransA, TransB, M, N, block_K, alpha, tmpA, lda, tmpB, ldb, Ires, ldc);
      }
      for(i = 0; i < M; i++){
        for(j = 0; j < N; j++){
          switch(Order){
            case 'r':
            case 'R':
              idxd_zziconv_sub(fold, Ires + (i * ldc + j) * idxd_zinum(fold), res + i * ldc + j);
              break;
            default:
              idxd_zziconv_sub(fold, Ires + (j * ldc + i) * idxd_zinum(fold), res + j * ldc + i);
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
              printf("reproBLAS_rzgemm(A, X, Y)[num_blocks=%d] = %g + %gi != %g + %gi\n", num_blocks, creal(res[i * ldc + j]), cimag(res[i * ldc + j]), creal(ref[i * ldc + j]), cimag(ref[i * ldc + j]));
              return 1;
            }
            break;
          default:
            if(res[j * ldc + i] != ref[j * ldc + i]){
              printf("reproBLAS_rzgemm(A, X, Y)[num_blocks=%d] = %g + %gi != %g + %gi\n", num_blocks, creal(res[j * ldc + i]), cimag(res[j * ldc + i]), creal(ref[j * ldc + i]), cimag(ref[j * ldc + i]));
              return 1;
            }
            break;
        }
      }
    }
    num_blocks *= 2;
  }
  free(res);
  free(Ires);
  return 0;
}

int matmat_fill_show_help(void){
  corroborate_rzgemm_options_initialize();

  opt_show_option(fold);
  opt_show_option(max_blocks);
  opt_show_option(shuffles);
  return 0;
}

const char* matmat_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  corroborate_rzgemm_options_initialize();

  opt_eval_option(argc, argv, &fold);
  opt_eval_option(argc, argv, &max_blocks);

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Corroborate rzgemm fold=%d", fold._int.value);
  return name_buffer;
}

int matmat_fill_test(int argc, char** argv, char Order, char TransA, char TransB, int M, int N, int K, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillB, double RealScaleB, double ImagScaleB, int ldb, double RealBeta, double ImagBeta, int FillC, double RealScaleC, double ImagScaleC, int ldc){
  (void)ImagAlpha;
  (void)ImagBeta;
  int rc = 0;
  int i;
  int j;

  corroborate_rzgemm_options_initialize();

  opt_eval_option(argc, argv, &fold);
  opt_eval_option(argc, argv, &max_blocks);
  opt_eval_option(argc, argv, &shuffles);

  util_random_seed();
  char NTransA;
  int opAM;
  int opAK;
  int opBK;
  int opBN;

  switch(TransA){
    case 'n':
    case 'N':
      opAM = M;
      opAK = K;
      NTransA = 't';
      break;
    default:
      opAM = K;
      opAK = M;
      NTransA = 'n';
      break;
  }

  switch(TransB){
    case 'n':
    case 'N':
      opBK = K;
      opBN = N;
      break;
    default:
      opBK = N;
      opBN = K;
      break;
  }


  double complex *A  = util_zmat_alloc(Order, opAM, opAK, lda);
  double complex *B  = util_zmat_alloc(Order, opBK, opBN, ldb);
  double complex *C  = util_zmat_alloc(Order, M, N, ldc);
  double complex betaC;
  double complex alpha = RealAlpha + I * ImagAlpha;
  double complex beta = RealBeta + I * ImagBeta;
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
  double_complex_indexed *CI = malloc(CNM * idxd_zisize(fold._int.value));

  int *P;

  util_zmat_fill(Order, NTransA, opAM, opAK, A, lda, FillA, RealScaleA, ImagScaleA);
  util_zmat_fill(Order, TransB, opBK, opBN, B, ldb, FillB, RealScaleB, ImagScaleB);
  util_zmat_fill(Order, 'n', M, N, C, ldc, FillC, RealScaleC, ImagScaleC);
  for(i = 0; i < M; i++){
    for(j = 0; j < N; j++){
      if(beta == 0.0){
        switch(Order){
          case 'r':
          case 'R':
            idxd_zisetzero(fold._int.value, CI + (i * ldc + j) * idxd_zinum(fold._int.value));
            break;
          default:
            idxd_zisetzero(fold._int.value, CI + (j * ldc + i) * idxd_zinum(fold._int.value));
            break;
        }
      }else if(beta == 1.0){
        switch(Order){
          case 'r':
          case 'R':
            idxd_zizconv(fold._int.value, C + i * ldc + j, CI + (i * ldc + j) * idxd_zinum(fold._int.value));
            break;
          default:
            idxd_zizconv(fold._int.value, C + j * ldc + i, CI + (j * ldc + i) * idxd_zinum(fold._int.value));
            break;
        }
      }else{
        switch(Order){
          case 'r':
          case 'R':
            betaC = C[i * ldc + j] * beta;
            idxd_zizconv(fold._int.value, &betaC, CI + (i * ldc + j) * idxd_zinum(fold._int.value));
            break;
          default:
            betaC = C[j * ldc + i] * beta;
            idxd_zizconv(fold._int.value, &betaC, CI + (j * ldc + i) * idxd_zinum(fold._int.value));
            break;
        }
      }
    }
  }
  double complex *ref  = (double complex*)malloc(CNM * sizeof(double complex));

  //compute with unpermuted data
  memcpy(ref, C, CNM * sizeof(double complex));

  wrap_ref_rzgemm(fold._int.value, Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, ref, ldc);

  rc = corroborate_rzgemm(fold._int.value, Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, C, CI, ldc, ref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(K);
  util_zmat_row_reverse(Order, NTransA, opAM, opAK, A, lda, P, 1);
  util_zmat_row_permute(Order, TransB, opBK, opBN, B, ldb, P, 1, NULL, 1);
  free(P);

  rc = corroborate_rzgemm(fold._int.value, Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, C, CI, ldc, ref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(K);
  util_zmat_row_sort(Order, NTransA, opAM, opAK, A, lda, P, 1, util_Increasing, 0);
  util_zmat_row_permute(Order, TransB, opBK, opBN, B, ldb, P, 1, NULL, 1);
  free(P);

  rc = corroborate_rzgemm(fold._int.value, Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, C, CI, ldc, ref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(K);
  util_zmat_row_sort(Order, NTransA, opAM, opAK, A, lda, P, 1, util_Decreasing, 0);
  util_zmat_row_permute(Order, TransB, opBK, opBN, B, ldb, P, 1, NULL, 1);
  free(P);

  rc = corroborate_rzgemm(fold._int.value, Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, C, CI, ldc, ref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(K);
  util_zmat_row_sort(Order, NTransA, opAM, opAK, A, lda, P, 1, util_Increasing_Magnitude, 0);
  util_zmat_row_permute(Order, TransB, opBK, opBN, B, ldb, P, 1, NULL, 1);
  free(P);

  rc = corroborate_rzgemm(fold._int.value, Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, C, CI, ldc, ref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(K);
  util_zmat_row_sort(Order, NTransA, opAM, opAK, A, lda, P, 1, util_Decreasing_Magnitude, 0);
  util_zmat_row_permute(Order, TransB, opBK, opBN, B, ldb, P, 1, NULL, 1);
  free(P);

  rc = corroborate_rzgemm(fold._int.value, Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, C, CI, ldc, ref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  for(i = 0; i < shuffles._int.value; i++){
    P = util_identity_permutation(K);
    util_zmat_row_shuffle(Order, NTransA, opAM, opAK, A, lda, P, 1);
    util_zmat_row_permute(Order, TransB, opBK, opBN, B, ldb, P, 1, NULL, 1);
    free(P);

    rc = corroborate_rzgemm(fold._int.value, Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, C, CI, ldc, ref, max_blocks._int.value);
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
