#include <indexedBLAS.h>
#include <indexed.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "rdblas1_wrapper.h"

#include "../common/test_vecvec_fill_header.h"

static opt_option func_type = {._named.header.type       = opt_named,
                               ._named.header.short_name = 'w',
                               ._named.header.long_name  = "w_type",
                               ._named.header.help       = "wrapped function type",
                               ._named.required          = 1,
                               ._named.n_names           = wrap_rdblas1_n_names,
                               ._named.names             = (char**)wrap_rdblas1_names,
                               ._named.descs             = (char**)wrap_rdblas1_descs,
                               ._named.value             = wrap_RDSUM};

int verify_rdblas1_reproducibility(int N, double* X, int incX, double* Y, int incY, int func, double ref, Idouble Iref, int max_num_blocks) {
  // GENERATE DATA
  int i;
  double res;
  Idouble Ires;
  int num_blocks = 1;

  int block_N = (N + num_blocks - 1) / num_blocks;

  num_blocks = 1;
  while (num_blocks < N && num_blocks <= max_num_blocks) {
    if (num_blocks == 1)
      res = (wrap_rdblas1_func(func))(N, X, incX, Y, incY);
    else {
      block_N =  (N + num_blocks - 1) / num_blocks;
      dISetZero(Ires);
      for (i = 0; i < N; i += block_N) {
        block_N = block_N < N - i ? block_N : (N-i);
        Idouble foo = (wrap_Idblas1_func(func))(block_N, X + i * incX, incX, Y + i * incY, incY);
        didiadd(&foo, &Ires, DEFAULT_FOLD);
      }
      res = ddiconv(&Ires, DEFAULT_FOLD);
    }
    if (res != ref) {
      printf("%s(X, Y)[num_blocks=%d,block_N=%d] = %g != %g\n", wrap_rdblas1_names[func], num_blocks, block_N, res, ref);
      if (num_blocks == 1) {
        Ires = (wrap_Idblas1_func(func))(N, X, incX, Y, incY);
      }
      printf("Ref I_double:\n");
      diprint(&Iref, DEFAULT_FOLD);
      printf("\nRes I_double:\n");
      diprint(&Ires, DEFAULT_FOLD);
      printf("\n");
      return 1;
    }
    num_blocks *= 2;
  }
  return 0;
}

int vecvec_fill_show_help(void){
  opt_show_option(func_type);
  return 0;
}

const char* vecvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  opt_eval_option(argc, argv, &func_type);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Verify %s reproducibility", wrap_rdblas1_names[func_type._named.value]);
  return name_buffer;
}

int vecvec_fill_test(int argc, char** argv, int N, int FillX, double ScaleX, double CondX, int incX, int FillY, double ScaleY, double CondY, int incY){
  int rc = 0;
  double ref;
  Idouble Iref;
  int max_num_blocks = 1024;

  util_random_seed();

  double *X = util_dvec_alloc(N, incX);
  double *Y = util_dvec_alloc(N, incY);
  int *P;

  opt_eval_option(argc, argv, &func_type);

  util_dvec_fill(N, X, incX, FillX, ScaleX, CondX);
  util_dvec_fill(N, Y, incY, FillY, ScaleY, CondY);

  //nrm2 doesn't make sense with more than 1 block.
  if(func_type._named.value == wrap_RDNRM2){
    max_num_blocks = 1;
  }

  //compute with unpermuted data
  ref  = (wrap_rdblas1_func(func_type._named.value))(N, X, incX, Y, incY);
  Iref = (wrap_Idblas1_func(func_type._named.value))(N, X, incX, Y, incY);

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

  free(X);
  free(Y);

  return rc;
}
