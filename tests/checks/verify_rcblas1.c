#include <indexedBLAS.h>
#include <idxd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "../../config.h"
#include "rcblas1_wrapper.h"

#include "../common/test_vecvec_fill_header.h"

static opt_option func_type;
static void verify_rcblas1_options_initialize(void){
  func_type._named.header.type       = opt_named;
  func_type._named.header.short_name = 'w';
  func_type._named.header.long_name  = "w_type";
  func_type._named.header.help       = "wrapped function type";
  func_type._named.required          = 1;
  func_type._named.n_names           = wrap_rcblas1_n_names;
  func_type._named.names             = (char**)wrap_rcblas1_names;
  func_type._named.descs             = (char**)wrap_rcblas1_descs;
  func_type._named.value             = wrap_RCSUM;
}

int verify_rcblas1_reproducibility(int N, float complex* X, int incX, float complex* Y, int incY, int func, float complex ref, float_complex_indexed *Iref, int max_num_blocks) {
  // GENERATE DATA
  int i;
  float complex res;
  float_complex_indexed *Ires = cialloc(SIDEFAULTFOLD);
  int num_blocks = 1;

  int block_N = (N + num_blocks - 1) / num_blocks;

  num_blocks = 1;
  while (num_blocks < N && num_blocks <= max_num_blocks) {
    if (num_blocks == 1)
      res = (wrap_rcblas1_func(func))(N, X, incX, Y, incY);
    else {
      block_N =  (N + num_blocks - 1) / num_blocks;
      cisetzero(SIDEFAULTFOLD, Ires);
      for (i = 0; i < N; i += block_N) {
        block_N = block_N < N - i ? block_N : (N-i);
        (wrap_ciblas1_func(func))(block_N, X + i * incX, incX, Y + i * incY, incY, Ires);
      }
      cciconv_sub(SIDEFAULTFOLD, Ires, &res);
    }
    if (res != ref) {
      printf("%s(X, Y)[num_blocks=%d,block_N=%d] = %g + %gi != %g + %gi\n", wrap_rcblas1_names[func], num_blocks, block_N, crealf(res), cimagf(res), crealf(ref), cimagf(ref));
      if (num_blocks == 1) {
        cisetzero(SIDEFAULTFOLD, Ires);
        (wrap_ciblas1_func(func))(N, X, incX, Y, incY, Ires);
      }
      printf("Ref I_float_Complex:\n");
      ciprint(SIDEFAULTFOLD, Iref);
      printf("\nRes I_float_Complex:\n");
      ciprint(SIDEFAULTFOLD, Ires);
      printf("\n");
      return 1;
    }
    num_blocks *= 2;
  }
  free(Ires);
  return 0;
}

int vecvec_fill_show_help(void){
  verify_rcblas1_options_initialize();

  opt_show_option(func_type);
  return 0;
}

const char* vecvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  verify_rcblas1_options_initialize();

  opt_eval_option(argc, argv, &func_type);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Verify %s reproducibility", wrap_rcblas1_names[func_type._named.value]);
  return name_buffer;
}

int vecvec_fill_test(int argc, char** argv, int N, int FillX, double RealScaleX, double ImagScaleX, int incX, int FillY, double RealScaleY, double ImagScaleY, int incY){
  int rc = 0;
  float complex ref;
  float_complex_indexed *Iref = cialloc(SIDEFAULTFOLD);
  int max_num_blocks = 1024;

  verify_rcblas1_options_initialize();

  util_random_seed();

  float complex *X = util_cvec_alloc(N, incX);
  float complex *Y = util_cvec_alloc(N, incY);
  int *P;

  opt_eval_option(argc, argv, &func_type);

  util_cvec_fill(N, X, incX, FillX, RealScaleX, ImagScaleX);
  util_cvec_fill(N, Y, incY, FillY, RealScaleY, ImagScaleY);

  //nrm2 doesn't make sense with more than 1 block.
  if(func_type._named.value == wrap_RSCNRM2){
    max_num_blocks = 1;
  }

  //compute with unpermuted data
  ref  = (wrap_rcblas1_func(func_type._named.value))(N, X, incX, Y, incY);
  cisetzero(SIDEFAULTFOLD, Iref);
  (wrap_ciblas1_func(func_type._named.value))(N, X, incX, Y, incY, Iref);

  P = util_identity_permutation(N);
  util_cvec_reverse(N, X, incX, P, 1);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_rcblas1_reproducibility(N, X, incX, Y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_sort(N, X, incX, P, 1, util_Increasing);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_rcblas1_reproducibility(N, X, incX, Y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_sort(N, X, incX, P, 1, util_Decreasing);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_rcblas1_reproducibility(N, X, incX, Y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_sort(N, X, incX, P, 1, util_Increasing_Magnitude);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_rcblas1_reproducibility(N, X, incX, Y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_sort(N, X, incX, P, 1, util_Decreasing_Magnitude);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_rcblas1_reproducibility(N, X, incX, Y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_shuffle(N, X, incX, P, 1);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_rcblas1_reproducibility(N, X, incX, Y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_shuffle(N, X, incX, P, 1);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_rcblas1_reproducibility(N, X, incX, Y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_shuffle(N, X, incX, P, 1);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_rcblas1_reproducibility(N, X, incX, Y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_shuffle(N, X, incX, P, 1);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_rcblas1_reproducibility(N, X, incX, Y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  free(Iref);
  free(X);
  free(Y);

  return rc;
}
