#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <indexed.h>
#include <indexedBLAS.h>
#include <reproBLAS.h>

#include "../common/test_opt.h"
#include "../common/test_matmat_fill_header.h"

#include "wrap_rdgemm.h"

static opt_option max_blocks;
static opt_option fold;

static void verify_rdgemm_options_initialize(void){
  max_blocks._int.header.type       = opt_int;
  max_blocks._int.header.short_name = 'B';
  max_blocks._int.header.long_name  = "blocks";
  max_blocks._int.header.help       = "maximum number of blocks";
  max_blocks._int.required          = 0;
  max_blocks._int.min               = 1;
  max_blocks._int.max               = INT_MAX;
  max_blocks._int.value             = 1024;

  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 2;
  fold._int.max               = DIMAXFOLD;
  fold._int.value             = DIDEFAULTFOLD;
}

int validate_internal_dgemm(int fold, char Order, char TransA, char TransB, int M, int N, int K, double alpha, double *A, int lda, double* B, int ldb, double beta, double *C, int ldc, double *ref) {

  int i;
  int j;
  int k;
  int num_blocks = 1;
  int block_K;

  double *res;
  double tmpres;
  double tmpref;

  switch(Order){
    case 'r':
    case 'R':
      res = malloc(M * ldc * sizeof(double));
      memcpy(res, C, M * ldc * sizeof(double));
      break;
    default:
      res = malloc(ldc * N * sizeof(double));
      memcpy(res, C, ldc * N * sizeof(double));
      break;
  }

  wrap_rdgemm(fold, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, res, ldc);
  for(i = 0; i < M; i++){
    for(j = 0; j < N; j++){
      switch(Order){
        case 'r':
        case 'R':
          tmpres = res[i * ldc + j];
          tmpref = ref[i * N + j];
          break;
        default:
          tmpres = res[j * ldc + i];
          tmpref = ref[j * M + i];
          break;
      }
      error = fabs(tmpres - tmpref);
      bound = wrap_rdgemm_bound(fold, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc, res, ref, i, j);
      if (!util_dsoftequals(tmpres, tmpref, bound)) {
        //TODO these error messages need to go to stderr for all tests.
        printf("rdgemm(A, X, Y) = %g != %g\n|%g - %g| = %g > %g\n", tmpres, tmpref, tmpres, tmpref, error, bound);
        return 1;
      }
    }
  }
  return 0;
}

int matmat_fill_show_help(void){
  verify_rdgemm_options_initialize();

  opt_show_option(fold);
  opt_show_option(max_blocks);
  opt_show_option(shuffles);
  return 0;
}

const char* matmat_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  verify_rdgemm_options_initialize();

  opt_eval_option(argc, argv, &fold);
  opt_eval_option(argc, argv, &max_blocks);

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Validate rdgemm internally fold=%d", fold._int.value);
  return name_buffer;
}

int matmat_fill_test(int argc, char** argv, char Order, char TransA, char TransB, int M, int N, int K, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillB, double RealScaleB, double ImagScaleB, int ldb, double RealBeta, double ImagBeta, int FillC, double RealScaleC, double ImagScaleC, int ldc){
  (void)ImagAlpha;
  (void)ImagBeta;
  int rc = 0;
  int i;
  int j;

  verify_rdgemm_options_initialize();

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
  int CN;
  switch(Order){
    case 'r':
    case 'R':
      CN = M * ldc;
      break;
    default:
      CN = ldc * N;
      break;
  }

  int *P;

  util_dmat_fill(Order, NTransA, M, K, A, lda, FillA, RealScaleA, ImagScaleA);
  util_dmat_fill(Order, TransB, K, N, B, ldb, FillB, RealScaleB, ImagScaleB);
  util_dmat_fill(Order, 'n', M, N, C, ldc, FillC, RealScaleC, ImagScaleC);
  double *ref  = (double*)wrap_rdgemm_result(Order, TransA, TransB, M, N, K, RealAlpha, ImagAlpha, FillA, RealScaleA, ImagScaleA, FillB, RealScaleB, ImagScaleB, RealBeta, ImagBeta, FillC, RealScaleC, ImagScaleC);

  P = util_identity_permutation(K);
  util_dmat_row_reverse(Order, NTransA, M, K, A, lda, P, 1);
  util_dmat_row_permute(Order, TransB, K, N, B, ldb, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_dgemm(fold._int.value, Order, TransA, TransB, M, N, K, RealAlpha, A, lda, B, ldb, RealBeta, C, ldc, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(K);
  util_dmat_row_sort(Order, NTransA, M, K, A, lda, P, 1, util_Increasing, 0);
  util_dmat_row_permute(Order, TransB, K, N, B, ldb, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_dgemm(fold._int.value, Order, TransA, TransB, M, N, K, RealAlpha, A, lda, B, ldb, RealBeta, C, ldc, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(K);
  util_dmat_row_sort(Order, NTransA, M, K, A, lda, P, 1, util_Decreasing, 0);
  util_dmat_row_permute(Order, TransB, K, N, B, ldb, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_dgemm(fold._int.value, Order, TransA, TransB, M, N, K, RealAlpha, A, lda, B, ldb, RealBeta, C, ldc, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(K);
  util_dmat_row_sort(Order, NTransA, M, K, A, lda, P, 1, util_Increasing_Magnitude, 0);
  util_dmat_row_permute(Order, TransB, K, N, B, ldb, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_dgemm(fold._int.value, Order, TransA, TransB, M, N, K, RealAlpha, A, lda, B, ldb, RealBeta, C, ldc, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(K);
  util_dmat_row_sort(Order, NTransA, M, K, A, lda, P, 1, util_Decreasing_Magnitude, 0);
  util_dmat_row_permute(Order, TransB, K, N, B, ldb, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_dgemm(fold._int.value, Order, TransA, TransB, M, N, K, RealAlpha, A, lda, B, ldb, RealBeta, C, ldc, ref);
  if(rc != 0){
    return rc;
  }

  free(A);
  free(B);
  free(C);
  free(ref);

  return rc;
}
