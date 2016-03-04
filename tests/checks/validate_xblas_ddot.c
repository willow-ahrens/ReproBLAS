#include <idxdBLAS.h>
#include <idxd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "../../config.h"
#include "wrap_daugsum.h"

#include "../common/test_vecvec_header.h"

static opt_option augsum_func;
static opt_option fold;
static opt_option norm;

static void validate_xblas_ddot_options_initialize(void){
  augsum_func._named.header.type       = opt_named;
  augsum_func._named.header.short_name = 'w';
  augsum_func._named.header.long_name  = "augsum_func";
  augsum_func._named.header.help       = "augmented summation function";
  augsum_func._named.required          = 0;
  augsum_func._named.n_names           = wrap_daugsum_func_n_names;
  augsum_func._named.names             = (char**)wrap_daugsum_func_names;
  augsum_func._named.descs             = (char**)wrap_daugsum_func_descs;
  augsum_func._named.value             = wrap_daugsum_RDDOT;

  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 2;
  fold._int.max               = idxd_DIMAXFOLD;
  fold._int.value             = DIDEFAULTFOLD;

  norm._int.header.type       = opt_int;
  norm._int.header.short_name = 'n';
  norm._int.header.long_name  = "norm";
  norm._int.header.help       = "-1 for close to underflow, 1 for close to overflow";
  norm._int.required          = 0;
  norm._int.min               = -1;
  norm._int.max               = 1;
  norm._int.value             = 0;
}

int validate_xblas_ddot(int fold, int N, double* X, int incX, double* Y, int incY, int func, double r, double ref) {
  // GENERATE DATA
  double res;
  double error;
  double bound;
  double_indexed *ires = idxd_dialloc(fold);

  idxd_didconv(fold, r, ires);
  (wrap_diaugsum_func(func))(fold, N, X, incX, Y, incY, ires);
  res = idxd_ddiconv(fold, ires);
  error = fabs(res - ref);
  bound = wrap_daugsum_bound(fold, N, func, X, incX, Y, incY, res, ref);
  if (!util_dsoftequals(res, ref, bound)) {
    //TODO these error messages need to go to stderr for all tests.
    printf("%s(X, Y) = %g != %g\n|%g - %g| = %g > %g\n", wrap_daugsum_func_names[func], res, ref, res, ref, error, bound);
    printf("\nres double_indexed:\n");
    idxd_diprint(fold, ires);
    printf("\n");
    return 1;
  }
  free(ires);
  return 0;
}

int vecvec_show_help(void){
  validate_xblas_ddot_options_initialize();

  opt_show_option(augsum_func);
  opt_show_option(fold);
  opt_show_option(norm);
  return 0;
}

const char* vecvec_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  validate_xblas_ddot_options_initialize();

  opt_eval_option(argc, argv, &augsum_func);
  opt_eval_option(argc, argv, &fold);
  opt_eval_option(argc, argv, &norm);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Validate %s internally fold=%d", wrap_daugsum_func_names[augsum_func._named.value], fold._int.value);
  return name_buffer;
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY){
  int rc = 0;
  double ref;
  double r;

  validate_xblas_ddot_options_initialize();

  util_random_seed();

  opt_eval_option(argc, argv, &augsum_func);
  opt_eval_option(argc, argv, &fold);
  opt_eval_option(argc, argv, &norm);

  double *X = util_dvec_alloc(N * 4, incX);
  double *Y = util_dvec_alloc(N * 4, incY);
  int *P;

  double *X_gen = util_dvec_alloc(N, 1);
  double *Y_gen = util_dvec_alloc(N, 1);
  double alpha = 1.0;
  double beta = 1.0;
  double ref_t;
  int seed = rand();
  util_xblas_ddot_fill(N, 0, 0, norm._int.value, 'n',
                  &alpha, 1, &beta, 1,
                  X_gen, Y_gen, &seed,
                  &r, &ref, &ref_t);
  util_dvec_dotsplit_twins(N, X_gen, 1, X, incX);
  util_dvec_dotsplit_couples(N, Y_gen, 1, Y, incY);
  free(X_gen);
  free(Y_gen);
  N *= 4;

  rc = validate_xblas_ddot(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, r, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_reverse(N, X, incX, P, 1);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_xblas_ddot(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, r, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_sort(N, X, incX, P, 1, util_Increasing);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_xblas_ddot(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, r, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_sort(N, X, incX, P, 1, util_Decreasing);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_xblas_ddot(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, r, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_sort(N, X, incX, P, 1, util_Increasing_Magnitude);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_xblas_ddot(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, r, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_sort(N, X, incX, P, 1, util_Decreasing_Magnitude);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_xblas_ddot(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, r, ref);
  if(rc != 0){
    return rc;
  }

  free(X);
  free(Y);

  return rc;
}
