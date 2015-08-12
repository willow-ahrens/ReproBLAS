#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <idxd.h>
#include <indexedBLAS.h>
#include <reproBLAS.h>

#include "../common/test_opt.h"
#include "../common/test_matvec_fill_header.h"

#include "wrap_rdgemv.h"

static opt_option fold;

static void validate_internal_rdgemv_options_initialize(void){
  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 2;
  fold._int.max               = DIMAXFOLD;
  fold._int.value             = DIDEFAULTFOLD;
}

int validate_internal_dgemv(int fold, char Order, char TransA, int M, int N, int NX, int NY, double alpha, double *A, int lda, double* X, int incX, double beta, double *Y, int incY, double *ref) {
  int i;
  double error;
  double bound;

  double *res = malloc(NY * incY * sizeof(double));

  memcpy(res, Y, NY * incY * sizeof(double));
  wrap_rdgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, beta, res, incY);
  for(i = 0; i < NY; i++){
    error = fabs(res[i * incY] - ref[i]);
    bound = wrap_rdgemv_bound(fold, Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY, res, ref, i);
    if (!util_dsoftequals(res[i * incY], ref[i], bound)) {
      //TODO these error messages need to go to stderr for all tests.
      printf("rdgemv(A, X, Y) = %g != %g\n|%g - %g| = %g > %g\n", res[i * incY], ref[i], res[i * incY], ref[i], error, bound);
      return 1;
    }
  }
  return 0;
}

int matvec_fill_show_help(void){
  validate_internal_rdgemv_options_initialize();

  opt_show_option(fold);
  return 0;
}

const char* matvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  validate_internal_rdgemv_options_initialize();

  opt_eval_option(argc, argv, &fold);

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Validate rdgemv internally fold=%d", fold._int.value);
  return name_buffer;
}

int matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillX, double RealScaleX, double ImagScaleX, int incX, double RealBeta, double ImagBeta, int FillY, double RealScaleY, double ImagScaleY, int incY){
  (void)ImagAlpha;
  (void)ImagBeta;
  int rc = 0;

  validate_internal_rdgemv_options_initialize();

  opt_eval_option(argc, argv, &fold);

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

  int *P;

  util_dmat_fill(Order, TransA, M, N, A, lda, FillA, RealScaleA, ImagScaleA);
  util_dvec_fill(NX, X, incX, FillX, RealScaleX, ImagScaleX);
  util_dvec_fill(NY, Y, incY, FillY, RealScaleY, ImagScaleY);
  double *ref = wrap_rdgemv_result(Order, TransA, M, N, RealAlpha, ImagAlpha, FillA, RealScaleA, ImagScaleA, A, lda, FillX, RealScaleX, ImagScaleX, X, incX, RealBeta, ImagBeta, Y, incY);

  rc = validate_internal_dgemv(fold._int.value, Order, TransA, M, N, NX, NY, RealAlpha, A, lda, X, incX, RealBeta, Y, incY, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_dvec_reverse(NX, X, incX, P, 1);
  util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_dgemv(fold._int.value, Order, TransA, M, N, NX, NY, RealAlpha, A, lda, X, incX, RealBeta, Y, incY, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_dvec_sort(NX, X, incX, P, 1, util_Increasing);
  util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_dgemv(fold._int.value, Order, TransA, M, N, NX, NY, RealAlpha, A, lda, X, incX, RealBeta, Y, incY, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_dvec_sort(NX, X, incX, P, 1, util_Decreasing);
  util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_dgemv(fold._int.value, Order, TransA, M, N, NX, NY, RealAlpha, A, lda, X, incX, RealBeta, Y, incY, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_dvec_sort(NX, X, incX, P, 1, util_Increasing_Magnitude);
  util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_dgemv(fold._int.value, Order, TransA, M, N, NX, NY, RealAlpha, A, lda, X, incX, RealBeta, Y, incY, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_dvec_sort(NX, X, incX, P, 1, util_Decreasing_Magnitude);
  util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_dgemv(fold._int.value, Order, TransA, M, N, NX, NY, RealAlpha, A, lda, X, incX, RealBeta, Y, incY, ref);
  if(rc != 0){
    return rc;
  }

  free(A);
  free(X);
  free(Y);
  free(ref);

  return rc;
}
