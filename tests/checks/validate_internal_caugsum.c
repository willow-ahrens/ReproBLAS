#include <indexedBLAS.h>
#include <indexed.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "../../config.h"
#include "wrap_caugsum.h"

#include "../common/test_vecvec_fill_header.h"

static opt_option augsum_func;
static opt_option fold;

static void validate_internal_caugsum_options_initialize(void){
  augsum_func._named.header.type       = opt_named;
  augsum_func._named.header.short_name = 'w';
  augsum_func._named.header.long_name  = "augsum_func";
  augsum_func._named.header.help       = "augmented summation function";
  augsum_func._named.required          = 1;
  augsum_func._named.n_names           = wrap_caugsum_func_n_names;
  augsum_func._named.names             = (char**)wrap_caugsum_func_names;
  augsum_func._named.descs             = (char**)wrap_caugsum_func_descs;
  augsum_func._named.value             = wrap_caugsum_RCSUM;

  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 2;
  fold._int.max               = SIMAXFOLD;
  fold._int.value             = SIDEFAULTFOLD;
}

int validate_internal_caugsum(int fold, int N, float complex* X, int incX, float complex* Y, int incY, int func, float complex ref) {
  // GENERATE DATA
  float complex res;
  float complex error;
  float complex bound;
  float_complex_indexed *ires = cialloc(fold);

  res = (wrap_caugsum_func(func))(fold, N, X, incX, Y, incY);
  error = res - ref;
  bound = wrap_caugsum_bound(fold, N, func, X, incX, Y, incY, res, ref);
  if (!util_csoftequals(res, ref, bound)) {
    //TODO these error messages need to go to stderr for all tests.
    printf("%s(X, Y) = %g + %gi != %g + %gi\n|%g - %g| = %g > %g and/or |%gi - %gi| = %g > %g\n", wrap_caugsum_func_names[func], crealf(res), cimagf(res), crealf(ref), cimagf(ref), crealf(res), crealf(ref), fabsf(crealf(error)), crealf(bound), cimagf(res), cimagf(ref), fabsf(cimagf(error)), cimagf(bound));
    cisetzero(fold, ires);
    (wrap_ciaugsum_func(func))(fold, N, X, incX, Y, incY, ires);
    printf("\nres float_complex_indexed:\n");
    ciprint(fold, ires);
    printf("\n");
    return 1;
  }
  free(ires);
  return 0;
}

int vecvec_fill_show_help(void){
  validate_internal_caugsum_options_initialize();

  opt_show_option(augsum_func);
  opt_show_option(fold);
  return 0;
}

const char* vecvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  validate_internal_caugsum_options_initialize();

  opt_eval_option(argc, argv, &augsum_func);
  opt_eval_option(argc, argv, &fold);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Validate %s internally fold=%d", wrap_caugsum_func_names[augsum_func._named.value], fold._int.value);
  return name_buffer;
}

int vecvec_fill_test(int argc, char** argv, int N, int FillX, double RealScaleX, double ImagScaleX, int incX, int FillY, double RealScaleY, double ImagScaleY, int incY){
  int rc = 0;
  float complex ref;

  validate_internal_caugsum_options_initialize();

  util_random_seed();

  opt_eval_option(argc, argv, &augsum_func);
  opt_eval_option(argc, argv, &fold);

  float complex *X = util_cvec_alloc(N, incX);
  float complex *Y = util_cvec_alloc(N, incY);
  int *P;

  util_cvec_fill(N, X, incX, FillX, RealScaleX, ImagScaleX);
  util_cvec_fill(N, Y, incY, FillY, RealScaleY, ImagScaleY);

  ref = wrap_caugsum_result(N, augsum_func._named.value, FillX, RealScaleX, ImagScaleX, FillY, RealScaleY, ImagScaleY);

  rc = validate_internal_caugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_reverse(N, X, incX, P, 1);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_caugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_sort(N, X, incX, P, 1, util_Increasing);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_caugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_sort(N, X, incX, P, 1, util_Decreasing);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_caugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_sort(N, X, incX, P, 1, util_Increasing_Magnitude);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_caugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_sort(N, X, incX, P, 1, util_Decreasing_Magnitude);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_caugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_shuffle(N, X, incX, P, 1);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_caugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_shuffle(N, X, incX, P, 1);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_caugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_shuffle(N, X, incX, P, 1);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_caugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_shuffle(N, X, incX, P, 1);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_caugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  free(X);
  free(Y);

  return rc;
}
