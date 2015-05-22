#include <indexedBLAS.h>
#include <indexed.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "../../config.h"
#include "wrap_daugsum.h"

#include "../common/test_vecvec_fill_header.h"

static opt_option augsum_func;
static opt_option fold;

static void validate_internal_daugsum_options_initialize(void){
  augsum_func._named.header.type       = opt_named;
  augsum_func._named.header.short_name = 'w';
  augsum_func._named.header.long_name  = "augsum_func";
  augsum_func._named.header.help       = "augmented summation function";
  augsum_func._named.required          = 1;
  augsum_func._named.n_names           = wrap_daugsum_func_n_names;
  augsum_func._named.names             = (char**)wrap_daugsum_func_names;
  augsum_func._named.descs             = (char**)wrap_daugsum_func_descs;
  augsum_func._named.value             = wrap_daugsum_RDSUM;

  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 1;
  fold._int.max               = MAX_FOLD;
  fold._int.value             = DEFAULT_FOLD;
}

int validate_internal_daugsum(int fold, int N, double* X, int incX, double* Y, int incY, int func, double ref) {
  // GENERATE DATA
  double res;
  double_indexed *ires = dialloc(fold);

  res = (wrap_daugsum_func(func))(fold, N, X, incX, Y, incY);
  if (res != ref) {
    //TODO these error messages need to go to stderr for all tests.
    printf("%s(X, Y) = %g != %g\n", wrap_daugsum_func_names[func], res, ref);
    disetzero(fold, ires);
    (wrap_diaugsum_func(func))(fold, N, X, incX, Y, incY, ires);
    printf("\nres double_indexed:\n");
    diprint(fold, ires);
    printf("\n");
    return 1;
  }
  free(ires);
  return 0;
}

int vecvec_fill_show_help(void){
  validate_internal_daugsum_options_initialize();

  opt_show_option(augsum_func);
  opt_show_option(fold);
  return 0;
}

const char* vecvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  validate_internal_daugsum_options_initialize();

  opt_eval_option(argc, argv, &augsum_func);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Validate %s internally", wrap_daugsum_func_names[augsum_func._named.value]);
  return name_buffer;
}

int vecvec_fill_test(int argc, char** argv, int N, int FillX, double ScaleX, double CondX, int incX, int FillY, double ScaleY, double CondY, int incY){
  int rc = 0;
  double ref;

  validate_internal_daugsum_options_initialize();

  util_random_seed();

  opt_eval_option(argc, argv, &augsum_func);
  opt_eval_option(argc, argv, &fold);

  double *X = util_dvec_alloc(N, incX);
  double *Y = util_dvec_alloc(N, incY);
  int *P;

  util_dvec_fill(N, X, incX, FillX, ScaleX, CondX);
  util_dvec_fill(N, Y, incY, FillY, ScaleY, CondY);


  //compute with unpermuted data
  ref = wrap_daugsum_result(N, augsum_func._named.value, FillX, ScaleX, CondX, FillY, ScaleY, CondY);

  rc = validate_internal_daugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_reverse(N, X, incX, P, 1);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_daugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_sort(N, X, incX, P, 1, util_Increasing);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_daugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_sort(N, X, incX, P, 1, util_Decreasing);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_daugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_sort(N, X, incX, P, 1, util_Increasing_Magnitude);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_daugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_sort(N, X, incX, P, 1, util_Decreasing_Magnitude);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_daugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_shuffle(N, X, incX, P, 1);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_daugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_shuffle(N, X, incX, P, 1);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_daugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_shuffle(N, X, incX, P, 1);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_daugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_shuffle(N, X, incX, P, 1);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_internal_daugsum(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref);
  if(rc != 0){
    return rc;
  }

  free(X);
  free(Y);

  return rc;
}
