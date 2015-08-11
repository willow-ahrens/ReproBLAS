#include <indexedBLAS.h>
#include <indexed.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "../../config.h"
#include "wrap_zaugsum.h"

#include "../common/test_vecvec_fill_header.h"

static opt_option augsum_func;
static opt_option fold;

static void validate_internal_zaugsum_options_initialize(void){
  augsum_func._named.header.type       = opt_named;
  augsum_func._named.header.short_name = 'w';
  augsum_func._named.header.long_name  = "augsum_func";
  augsum_func._named.header.help       = "augmented summation function";
  augsum_func._named.required          = 1;
  augsum_func._named.n_names           = wrap_zaugsum_func_n_names;
  augsum_func._named.names             = (char**)wrap_zaugsum_func_names;
  augsum_func._named.descs             = (char**)wrap_zaugsum_func_descs;
  augsum_func._named.value             = wrap_zaugsum_RZSUM;

  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 2;
  fold._int.max               = DIMAXFOLD;
  fold._int.value             = DIDEFAULTFOLD;
}

int validate_internal_zaugsum(int fold, int N, double complex* X, int incX, double complex* Y, int incY, int func, double complex ref) {
  // GENERATE DATA
  double complex res;
  double complex error;
  double complex bound;
  double_complex_indexed *ires = zialloc(fold);

  res = (wrap_zaugsum_func(func))(fold, N, X, incX, Y, incY);
  error = res - ref;
  bound = wrap_zaugsum_bound(fold, N, func, X, incX, Y, incY, res, ref);
  if (!util_zsoftequals(res, ref, bound)) {
    //TODO these error messages need to go to stderr for all tests.
    printf("%s(X, Y) = %g + %gi != %g + %gi\n|%g - %g| = %g > %g and/or |%gi - %gi| = %g > %g\n", wrap_zaugsum_func_names[func], creal(res), cimag(res), creal(ref), cimag(ref), creal(res), creal(ref), fabs(creal(error)), creal(bound), cimag(res), cimag(ref), fabs(cimag(error)), cimag(bound));
    zisetzero(fold, ires);
    (wrap_ziaugsum_func(func))(fold, N, X, incX, Y, incY, ires);
    printf("\nres double_complex_indexed:\n");
    ziprint(fold, ires);
    printf("\n");
    return 1;
  }
  free(ires);
  return 0;
}

int vecvec_fill_show_help(void){
  validate_internal_zaugsum_options_initialize();

  opt_show_option(augsum_func);
  opt_show_option(fold);
  return 0;
}

const char* vecvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  validate_internal_zaugsum_options_initialize();

  opt_eval_option(argc, argv, &augsum_func);
  opt_eval_option(argc, argv, &fold);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Validate %s internally fold=%d", wrap_zaugsum_func_names[augsum_func._named.value], fold._int.value);
  return name_buffer;
}

int vecvec_fill_test(int argc, char** argv, int N, int FillX, double RealScaleX, double ImagScaleX, int incX, int FillY, double RealScaleY, double ImagScaleY, int incY){
  int rc = 0;
  double complex ref;

  validate_internal_zaugsum_options_initialize();

  util_random_seed();

  opt_eval_option(argc, argv, &augsum_func);
  opt_eval_option(argc, argv, &fold);

  double complex *X = util_zvec_alloc(N, incX);
  double complex *Y = util_zvec_alloc(N, incY);
  int *P;

  util_zvec_fill(N, X, incX, FillX, RealScaleX, ImagScaleX);
  util_zvec_fill(N, Y, incY, FillY, RealScaleY, ImagScaleY);

  ref = wrap_zaugsum_result(N, augsum_func._named.value, FillX, RealScaleX, ImagScaleX, FillY, RealScaleY, ImagScaleY);

  rc = validate_internal_zaugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_zvec_reverse(N, X, incX, P, 1);
  util_zvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_zaugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_zvec_sort(N, X, incX, P, 1, util_Increasing);
  util_zvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_zaugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_zvec_sort(N, X, incX, P, 1, util_Decreasing);
  util_zvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_zaugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_zvec_sort(N, X, incX, P, 1, util_Increasing_Magnitude);
  util_zvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_zaugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_zvec_sort(N, X, incX, P, 1, util_Decreasing_Magnitude);
  util_zvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_zaugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  free(X);
  free(Y);

  return rc;
}
