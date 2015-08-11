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
static opt_option max_blocks;
static opt_option shuffles;
static opt_option fold;

static void verify_daugsum_options_initialize(void){
  augsum_func._named.header.type       = opt_named;
  augsum_func._named.header.short_name = 'w';
  augsum_func._named.header.long_name  = "augsum_func";
  augsum_func._named.header.help       = "augmented summation function";
  augsum_func._named.required          = 1;
  augsum_func._named.n_names           = wrap_daugsum_func_n_names;
  augsum_func._named.names             = (char**)wrap_daugsum_func_names;
  augsum_func._named.descs             = (char**)wrap_daugsum_func_descs;
  augsum_func._named.value             = wrap_daugsum_RDSUM;

  max_blocks._int.header.type       = opt_int;
  max_blocks._int.header.short_name = 'B';
  max_blocks._int.header.long_name  = "blocks";
  max_blocks._int.header.help       = "maximum number of blocks";
  max_blocks._int.required          = 0;
  max_blocks._int.min               = 1;
  max_blocks._int.max               = INT_MAX;
  max_blocks._int.value             = 1024;

  shuffles._int.header.type       = opt_int;
  shuffles._int.header.short_name = 'S';
  shuffles._int.header.long_name  = "shuffles";
  shuffles._int.header.help       = "number of times to shuffle";
  shuffles._int.required          = 0;
  shuffles._int.min               = 0;
  shuffles._int.max               = INT_MAX;
  shuffles._int.value             = 5;

  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 2;
  fold._int.max               = DIMAXFOLD;
  fold._int.value             = DIDEFAULTFOLD;
}

int verify_daugsum_reproducibility(int fold, int N, double* X, int incX, double* Y, int incY, int func, double ref, double_indexed *iref, int max_num_blocks) {
  // GENERATE DATA
  int i;
  double res;
  double_indexed *ires = dialloc(fold);
  int num_blocks = 1;

  int block_N = (N + num_blocks - 1) / num_blocks;

  num_blocks = 1;
  while (num_blocks < N && num_blocks <= max_num_blocks) {
    if (num_blocks == 1)
      res = (wrap_daugsum_func(func))(fold, N, X, incX, Y, incY);
    else {
      block_N =  (N + num_blocks - 1) / num_blocks;
      disetzero(fold, ires);
      for (i = 0; i < N; i += block_N) {
        block_N = block_N < N - i ? block_N : (N-i);
        (wrap_diaugsum_func(func))(fold, block_N, X + i * incX, incX, Y + i * incY, incY, ires);
      }
      res = ddiconv(fold, ires);
    }
    if (res != ref) {
      printf("%s(X, Y)[num_blocks=%d,block_N=%d] = %g != %g\n", wrap_daugsum_func_names[func], num_blocks, block_N, res, ref);
      if (num_blocks == 1) {
        disetzero(fold, ires);
        (wrap_diaugsum_func(func))(fold, N, X, incX, Y, incY, ires);
      }
      printf("ref double_indexed:\n");
      diprint(fold, iref);
      printf("\nres double_indexed:\n");
      diprint(fold, ires);
      printf("\n");
      return 1;
    }
    num_blocks *= 2;
  }
  free(ires);
  return 0;
}

int vecvec_fill_show_help(void){
  verify_daugsum_options_initialize();

  opt_show_option(augsum_func);
  opt_show_option(max_blocks);
  opt_show_option(shuffles);
  opt_show_option(fold);
  return 0;
}

const char* vecvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  verify_daugsum_options_initialize();

  opt_eval_option(argc, argv, &augsum_func);
  opt_eval_option(argc, argv, &fold);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Verify %s reproducibility fold=%d", wrap_daugsum_func_names[augsum_func._named.value], fold._int.value);
  return name_buffer;
}

int vecvec_fill_test(int argc, char** argv, int N, int FillX, double RealScaleX, double ImagScaleX, int incX, int FillY, double RealScaleY, double ImagScaleY, int incY){
  int rc = 0;
  int i;
  double ref;
  double_indexed *iref;
  int max_num_blocks;

  verify_daugsum_options_initialize();

  util_random_seed();

  opt_eval_option(argc, argv, &augsum_func);
  opt_eval_option(argc, argv, &max_blocks);
  opt_eval_option(argc, argv, &shuffles);
  opt_eval_option(argc, argv, &fold);

  iref = dialloc(fold._int.value);

  double *X = util_dvec_alloc(N, incX);
  double *Y = util_dvec_alloc(N, incY);
  int *P;

  max_num_blocks = max_blocks._int.value;
  //nrm2 doesn't make sense with more than 1 block.
  if(augsum_func._named.value == wrap_daugsum_RDNRM2){
    max_num_blocks = 1;
  }

  util_dvec_fill(N, X, incX, FillX, RealScaleX, ImagScaleX);
  util_dvec_fill(N, Y, incY, FillY, RealScaleY, ImagScaleY);


  //compute with unpermuted data
  ref  = (wrap_daugsum_func(augsum_func._named.value))(fold._int.value, N, X, incX, Y, incY);
  disetzero(fold._int.value, iref);
  (wrap_diaugsum_func(augsum_func._named.value))(fold._int.value, N, X, incX, Y, incY, iref);

  P = util_identity_permutation(N);
  util_dvec_reverse(N, X, incX, P, 1);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_daugsum_reproducibility(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref, iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_sort(N, X, incX, P, 1, util_Increasing);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_daugsum_reproducibility(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref, iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_sort(N, X, incX, P, 1, util_Decreasing);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_daugsum_reproducibility(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref, iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_sort(N, X, incX, P, 1, util_Increasing_Magnitude);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_daugsum_reproducibility(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref, iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_dvec_sort(N, X, incX, P, 1, util_Decreasing_Magnitude);
  util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_daugsum_reproducibility(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref, iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  for(i = 0; i < shuffles._int.value; i++){
    P = util_identity_permutation(N);
    util_dvec_shuffle(N, X, incX, P, 1);
    util_dvec_permute(N, Y, incY, P, 1, NULL, 1);
    free(P);

    rc = verify_daugsum_reproducibility(fold._int.value, N, X, incX, Y, incY, augsum_func._named.value, ref, iref, max_num_blocks);
    if(rc != 0){
      return rc;
    }
  }

  free(iref);
  free(X);
  free(Y);

  return rc;
}
