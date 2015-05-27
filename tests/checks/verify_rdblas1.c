#include <indexedBLAS.h>
#include <indexed.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "rdblas1_wrapper.h"

#include "../common/test_vecvec_fill_header.h"

static opt_option func_type;
static void verify_rdblas1_options_initialize(void){
  func_type._named.header.type       = opt_named;
  func_type._named.header.short_name = 'w';
  func_type._named.header.long_name  = "w_type";
  func_type._named.header.help       = "wrapped function type";
  func_type._named.required          = 1;
  func_type._named.n_names           = wrap_rdblas1_n_names;
  func_type._named.names             = (char**)wrap_rdblas1_names;
  func_type._named.descs             = (char**)wrap_rdblas1_descs;
  func_type._named.value             = wrap_RDSUM;
}

int verify_rdblas1_reproducibility(int N, double* X, int incX, double* Y, int incY, int func, double ref, double_indexed *Iref, int max_num_blocks) {
  // GENERATE DATA
  int i;
  double res;
  double_indexed *Ires = dialloc(DEFAULT_FOLD);
  int num_blocks = 1;

  int block_N = (N + num_blocks - 1) / num_blocks;

  num_blocks = 1;
  while (num_blocks < N && num_blocks <= max_num_blocks) {
    if (num_blocks == 1)
      res = (wrap_rdblas1_func(func))(N, X, incX, Y, incY);
    else {
      block_N =  (N + num_blocks - 1) / num_blocks;
      disetzero(DEFAULT_FOLD, Ires);
      for (i = 0; i < N; i += block_N) {
        block_N = block_N < N - i ? block_N : (N-i);
        (wrap_diblas1_func(func))(block_N, X + i * incX, incX, Y + i * incY, incY, Ires);
      }
      res = ddiconv(DEFAULT_FOLD, Ires);
    }
    if (res != ref) {
      printf("%s(X, Y)[num_blocks=%d,block_N=%d] = %g != %g\n", wrap_rdblas1_names[func], num_blocks, block_N, res, ref);
      if (num_blocks == 1) {
        disetzero(DEFAULT_FOLD, Ires);
        (wrap_diblas1_func(func))(N, X, incX, Y, incY, Ires);
      }
      printf("Ref I_double:\n");
      diprint(DEFAULT_FOLD, Iref);
      printf("\nRes I_double:\n");
      diprint(DEFAULT_FOLD, Ires);
      printf("\n");
      return 1;
    }
    num_blocks *= 2;
  }
  free(Ires);
  return 0;
}

int vecvec_fill_show_help(void){
  verify_rdblas1_options_initialize();

  opt_show_option(func_type);
  return 0;
}

const char* vecvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  verify_rdblas1_options_initialize();

  opt_eval_option(argc, argv, &func_type);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Verify %s reproducibility", wrap_rdblas1_names[func_type._named.value]);
  return name_buffer;
}

int vecvec_fill_test(int argc, char** argv, int N, int FillX, double RealScaleX, double ImagScaleX, int incX, int FillY, double RealScaleY, double ImagScaleY, int incY){
  int rc = 0;
  double ref;
  double_indexed *Iref = dialloc(DEFAULT_FOLD);
  int max_num_blocks = 1024;

  verify_rdblas1_options_initialize();

  util_random_seed();

  double *X = util_dvec_alloc(N, incX);
  double *Y = util_dvec_alloc(N, incY);
  int *P;

  opt_eval_option(argc, argv, &func_type);

  util_dvec_fill(N, X, incX, FillX, RealScaleX, ImagScaleX);
  util_dvec_fill(N, Y, incY, FillY, RealScaleY, ImagScaleY);

  //nrm2 doesn't make sense with more than 1 block.
  if(func_type._named.value == wrap_RDNRM2){
    max_num_blocks = 1;
  }

  //compute with unpermuted data
  ref  = (wrap_rdblas1_func(func_type._named.value))(N, X, incX, Y, incY);
  disetzero(DEFAULT_FOLD, Iref);
  (wrap_diblas1_func(func_type._named.value))(N, X, incX, Y, incY, Iref);

  P = util_identity_permutation(N);
  util_dvec_reverse(N, X, incX, P, 1);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_rdblas1_reproducibility(N, X, incX, Y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_sort(N, X, incX, P, 1, util_Increasing);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_rdblas1_reproducibility(N, X, incX, Y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_sort(N, X, incX, P, 1, util_Decreasing);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_rdblas1_reproducibility(N, X, incX, Y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_sort(N, X, incX, P, 1, util_Increasing_Magnitude);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_rdblas1_reproducibility(N, X, incX, Y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_sort(N, X, incX, P, 1, util_Decreasing_Magnitude);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_rdblas1_reproducibility(N, X, incX, Y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_shuffle(N, X, incX, P, 1);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_rdblas1_reproducibility(N, X, incX, Y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_shuffle(N, X, incX, P, 1);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_rdblas1_reproducibility(N, X, incX, Y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_shuffle(N, X, incX, P, 1);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_rdblas1_reproducibility(N, X, incX, Y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_shuffle(N, X, incX, P, 1);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_rdblas1_reproducibility(N, X, incX, Y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  free(Iref);
  free(X);
  free(Y);

  return rc;
}
