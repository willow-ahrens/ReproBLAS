#include <indexedBLAS.h>
#include <indexed.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "../common/test_opt.h"
#include "../../config.h"

#include "../common/test_vecvec_fill_header.h"

static opt_option max_blocks;
static opt_option fold;

static void verify_sisssq_options_initialize(void){
  max_blocks._int.header.type       = opt_int;
  max_blocks._int.header.short_name = 'B';
  max_blocks._int.header.long_name  = "blocks";
  max_blocks._int.header.help       = "maximum number of blocks";
  max_blocks._int.required          = 0;
  max_blocks._int.min               = 1;
  max_blocks._int.max               = INT_MAX;
  max_blocks._int.value             = 1024;

  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 2;
  fold._int.max               = MAX_FOLD;
  fold._int.value             = DEFAULT_FOLD;
}

int verify_sisssq_reproducibility(int fold, int N, float* X, int incX, float refscl, float refssq, float_indexed *iref, int max_num_blocks) {
  // GENERATE DATA
  int i;
  float resscl = 0.0;
  float resssq;
  float_indexed *ires = sialloc(fold);
  int num_blocks = 1;

  int block_N = (N + num_blocks - 1) / num_blocks;

  num_blocks = 1;
  while (num_blocks < N && num_blocks <= max_num_blocks) {
    sisetzero(fold, ires);
    if (num_blocks == 1)
      resscl = sisssq(fold, N, X, incX, resscl, ires);
    else {
      block_N = (N + num_blocks - 1) / num_blocks;
      for (i = 0; i < N; i += block_N) {
        block_N = block_N < N - i ? block_N : (N-i);
        resscl = sisssq(fold, block_N, X + i * incX, incX, resscl, ires);
      }
    }
    resssq = ssiconv(fold, ires);
    if (resssq != refssq || resscl != refscl) {
      printf("sisssq(X)[num_blocks=%d,block_N=%d] = %g * (%g) != %g * (%g)\n", num_blocks, block_N, resscl, resssq, refscl, refssq);
      printf("ref float_indexed:\n");
      siprint(fold, iref);
      printf("\nres float_indexed:\n");
      siprint(fold, ires);
      printf("\n");
      return 1;
    }
    num_blocks *= 2;
  }
  free(ires);
  return 0;
}

int vecvec_fill_show_help(void){
  verify_sisssq_options_initialize();

  opt_show_option(max_blocks);
  opt_show_option(fold);
  return 0;
}

const char* vecvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  verify_sisssq_options_initialize();

  opt_eval_option(argc, argv, &fold);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Verify sisssq reproducibility fold=%d", fold._int.value);
  return name_buffer;
}

int vecvec_fill_test(int argc, char** argv, int N, int FillX, double RealScaleX, double ImagScaleX, int incX, int FillY, double RealScaleY, double ImagScaleY, int incY){
  int rc = 0;
  float refscl = 0.0;
  float refssq;
  float_indexed *iref;
  int max_num_blocks;

  verify_sisssq_options_initialize();

  util_random_seed();

  opt_eval_option(argc, argv, &max_blocks);
  opt_eval_option(argc, argv, &fold);

  iref = sialloc(fold._int.value);

  float *X = util_svec_alloc(N, incX);
  float *Y = util_svec_alloc(N, incY);
  int *P;


  max_num_blocks = max_blocks._int.value;

  util_svec_fill(N, X, incX, FillX, RealScaleX, ImagScaleX);
  util_svec_fill(N, Y, incY, FillY, RealScaleY, ImagScaleY);


  //compute with unpermuted data
  sisetzero(fold._int.value, iref);
  refscl = sisssq(fold._int.value, N, X, incX, refscl, iref);
  refssq = ssiconv(fold._int.value, iref);

  P = util_identity_permutation(N);
  util_svec_reverse(N, X, incX, P, 1);
  util_svec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_sisssq_reproducibility(fold._int.value, N, X, incX, refscl, refssq, iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_svec_sort(N, X, incX, P, 1, util_Increasing);
  util_svec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_sisssq_reproducibility(fold._int.value, N, X, incX, refscl, refssq, iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_svec_sort(N, X, incX, P, 1, util_Decreasing);
  util_svec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_sisssq_reproducibility(fold._int.value, N, X, incX, refscl, refssq, iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_svec_sort(N, X, incX, P, 1, util_Increasing_Magnitude);
  util_svec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_sisssq_reproducibility(fold._int.value, N, X, incX, refscl, refssq, iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_svec_sort(N, X, incX, P, 1, util_Decreasing_Magnitude);
  util_svec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_sisssq_reproducibility(fold._int.value, N, X, incX, refscl, refssq, iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_svec_shuffle(N, X, incX, P, 1);
  util_svec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_sisssq_reproducibility(fold._int.value, N, X, incX, refscl, refssq, iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_svec_shuffle(N, X, incX, P, 1);
  util_svec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_sisssq_reproducibility(fold._int.value, N, X, incX, refscl, refssq, iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_svec_shuffle(N, X, incX, P, 1);
  util_svec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_sisssq_reproducibility(fold._int.value, N, X, incX, refscl, refssq, iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_svec_shuffle(N, X, incX, P, 1);
  util_svec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_sisssq_reproducibility(fold._int.value, N, X, incX, refscl, refssq, iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  free(iref);
  free(X);
  free(Y);

  return rc;
}
