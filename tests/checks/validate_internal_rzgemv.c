#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <indexed.h>
#include <indexedBLAS.h>
#include <reproBLAS.h>

#include "../common/test_opt.h"
#include "../common/test_matvec_fill_header.h"

#include "wrap_rzgemv.h"

static opt_option fold;

static void validate_internal_rzgemv_options_initialize(void){
  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 2;
  fold._int.max               = DIMAXFOLD;
  fold._int.value             = DIDEFAULTFOLD;
}

int validate_internal_zgemv(int fold, char Order, char TransA, int M, int N, int NX, int NY, double complex *alpha, double complex *A, int lda, double complex *X, int incX, double complex *beta, double complex *Y, int incY, double complex *ref) {
  int i;
  double complex error;
  double complex bound;

  double complex *res = malloc(NY * incY * sizeof(double complex));

  memcpy(res, Y, NY * incY * sizeof(double complex));
  wrap_rzgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, beta, res, incY);
  for(i = 0; i < NY; i++){
    error = res[i * incY] - ref[i];
    bound = wrap_rzgemv_bound(fold, Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY, res, ref, i);
    if (!util_zsoftequals(res[i * incY], ref[i], bound)) {
      //TODO these error messages need to go to stderr for all tests.
      printf("rzgemv(A, X, Y) = %g + %gi != %g + %gi\n|%g - %g| = %g > %g and/or |%gi - %gi| = %g > %g\n", creal(res[i * incY]), cimag(res[i * incY]), creal(ref[i]), cimag(ref[i]), creal(res[i * incY]), creal(ref[i]), fabs(creal(error)), creal(bound), cimag(res[i * incY]), cimag(ref[i]), fabs(cimag(error)), cimag(bound));
      return 1;
    }
  }
  return 0;
}

int matvec_fill_show_help(void){
  validate_internal_rzgemv_options_initialize();

  opt_show_option(fold);
  return 0;
}

const char* matvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  validate_internal_rzgemv_options_initialize();

  opt_eval_option(argc, argv, &fold);

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Validate rzgemv internally fold=%d", fold._int.value);
  return name_buffer;
}

int matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillX, double RealScaleX, double ImagScaleX, int incX, double RealBeta, double ImagBeta, int FillY, double RealScaleY, double ImagScaleY, int incY){
  int rc = 0;

  validate_internal_rzgemv_options_initialize();

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

  double complex *A  = util_zmat_alloc(Order, M, N, lda);
  double complex *X  = util_zvec_alloc(NX, incX);
  double complex *Y  = util_zvec_alloc(NY, incY);
  double complex alpha = RealAlpha + I * ImagAlpha;
  double complex beta = RealBeta + I * ImagBeta;

  int *P;

  util_zmat_fill(Order, TransA, M, N, A, lda, FillA, RealScaleA, ImagScaleA);
  util_zvec_fill(NX, X, incX, FillX, RealScaleX, ImagScaleX);
  util_zvec_fill(NY, Y, incY, FillY, RealScaleY, ImagScaleY);
  double complex *ref = wrap_rzgemv_result(Order, TransA, M, N, RealAlpha, ImagAlpha, FillA, RealScaleA, ImagScaleA, A, lda, FillX, RealScaleX, ImagScaleX, X, incX, RealBeta, ImagBeta, Y, incY);

  rc = validate_internal_zgemv(fold._int.value, Order, TransA, M, N, NX, NY, &alpha, A, lda, X, incX, &beta, Y, incY, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_zvec_reverse(NX, X, incX, P, 1);
  util_zmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_zgemv(fold._int.value, Order, TransA, M, N, NX, NY, &alpha, A, lda, X, incX, &beta, Y, incY, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_zvec_sort(NX, X, incX, P, 1, util_Increasing);
  util_zmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_zgemv(fold._int.value, Order, TransA, M, N, NX, NY, &alpha, A, lda, X, incX, &beta, Y, incY, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_zvec_sort(NX, X, incX, P, 1, util_Decreasing);
  util_zmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_zgemv(fold._int.value, Order, TransA, M, N, NX, NY, &alpha, A, lda, X, incX, &beta, Y, incY, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_zvec_sort(NX, X, incX, P, 1, util_Increasing_Magnitude);
  util_zmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_zgemv(fold._int.value, Order, TransA, M, N, NX, NY, &alpha, A, lda, X, incX, &beta, Y, incY, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_zvec_sort(NX, X, incX, P, 1, util_Decreasing_Magnitude);
  util_zmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_zgemv(fold._int.value, Order, TransA, M, N, NX, NY, &alpha, A, lda, X, incX, &beta, Y, incY, ref);
  if(rc != 0){
    return rc;
  }

  free(A);
  free(X);
  free(Y);
  free(ref);

  return rc;
}
