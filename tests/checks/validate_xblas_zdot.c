#include <idxdBLAS.h>
#include <idxd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "../../config.h"
#include "wrap_zaugsum.h"

#include "../common/test_vecvec_header.h"

static opt_option augsum_func;
static opt_option fold;
static opt_option norm;

static void validate_xblas_zdot_options_initialize(void){
  augsum_func._named.header.type       = opt_named;
  augsum_func._named.header.short_name = 'w';
  augsum_func._named.header.long_name  = "augsum_func";
  augsum_func._named.header.help       = "augmented summation function";
  augsum_func._named.required          = 0;
  augsum_func._named.n_names           = wrap_zaugsum_func_n_names;
  augsum_func._named.names             = (char**)wrap_zaugsum_func_names;
  augsum_func._named.descs             = (char**)wrap_zaugsum_func_descs;
  augsum_func._named.value             = wrap_zaugsum_RZDOTU;

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

int validate_xblas_zdot(int fold, int N, double complex* X, int incX, double complex* Y, int incY, int func, double complex r, double complex ref) {
  // GENERATE DATA
  double complex res;
  double complex error;
  double complex bound;
  double_complex_indexed *ires = idxd_zialloc(fold);

  idxd_zizconv(fold, (double*)&r, ires);
  (wrap_ziaugsum_func(func))(fold, N, X, incX, Y, incY, ires);
  idxd_zziconv_sub(fold, ires, (double*)&res);
  error = res - ref;
  bound = wrap_zaugsum_bound(fold, N, func, X, incX, Y, incY, res, ref);
  if (!util_zsoftequals(res, ref, bound)) {
    //TODO these error messages need to go to stderr for all tests.
    printf("%s(X, Y) = %g + %gi != %g + %gi\n|%g - %g| = %g > %g and/or |%gi - %gi| = %g > %g\n", wrap_zaugsum_func_names[func], creal(res), cimag(res), creal(ref), cimag(ref), creal(res), creal(ref), fabs(creal(error)), creal(bound), cimag(res), cimag(ref), fabs(cimag(error)), cimag(bound));
    printf("\nres double_complex_indexed:\n");
    idxd_ziprint(fold, ires);
    printf("\n");
    return 1;
  }
  free(ires);
  return 0;
}

int vecvec_show_help(void){
  validate_xblas_zdot_options_initialize();

  opt_show_option(augsum_func);
  opt_show_option(fold);
  opt_show_option(norm);
  return 0;
}

const char* vecvec_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  validate_xblas_zdot_options_initialize();

  opt_eval_option(argc, argv, &augsum_func);
  opt_eval_option(argc, argv, &fold);
  opt_eval_option(argc, argv, &norm);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Validate %s internally fold=%d", wrap_zaugsum_func_names[augsum_func._named.value], fold._int.value);
  return name_buffer;
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY){
  int rc = 0;
  double complex ref;
  double complex r;

  validate_xblas_zdot_options_initialize();

  util_random_seed();

  opt_eval_option(argc, argv, &augsum_func);
  opt_eval_option(argc, argv, &fold);
  opt_eval_option(argc, argv, &norm);

  double complex *X = util_zvec_alloc(N, incX * 4);
  double complex *Y = util_zvec_alloc(N, incY * 4);
  int *P;

  double complex *X_gen = util_zvec_alloc(N, 1);
  double complex *Y_gen = util_zvec_alloc(N, 1);
  double complex alpha = 1.0;
  double complex beta = 1.0;
  double complex ref_t;
  int seed = rand();
  char conj = 'n';
  if(augsum_func._named.value == wrap_zaugsum_RZDOTC){
    conj = 'c';
  }
  util_xblas_zdot_fill(N, 0, 0, norm._int.value, conj,
                  &alpha, 1, &beta, 1,
                  X_gen, Y_gen, &seed,
                  &r, (double*)&ref, (double*)&ref_t);
  util_zvec_dotsplit_twins(N, X_gen, 1, X, incX);
  util_zvec_dotsplit_couples(N, Y_gen, 1, Y, incY);
  free(X_gen);
  free(Y_gen);
  N *= 4;

  rc = validate_xblas_zdot(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, r, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_zvec_reverse(N, X, incX, P, 1);
  util_zvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_xblas_zdot(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, r, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_zvec_sort(N, X, incX, P, 1, util_Increasing);
  util_zvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_xblas_zdot(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, r, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_zvec_sort(N, X, incX, P, 1, util_Decreasing);
  util_zvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_xblas_zdot(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, r, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_zvec_sort(N, X, incX, P, 1, util_Increasing_Magnitude);
  util_zvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_xblas_zdot(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, r, ref);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_zvec_sort(N, X, incX, P, 1, util_Decreasing_Magnitude);
  util_zvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = validate_xblas_zdot(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, r, ref);
  if(rc != 0){
    return rc;
  }

  free(X);
  free(Y);

  return rc;
}
