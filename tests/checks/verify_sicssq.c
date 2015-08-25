#include <idxdBLAS.h>
#include <idxd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "../common/test_opt.h"
#include "../../config.h"

#include "../common/test_vecvec_fill_header.h"

static opt_option max_blocks;
static opt_option shuffles;
static opt_option fold;

static void verify_sicssq_options_initialize(void){
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
  fold._int.max               = SIMAXFOLD;
  fold._int.value             = SIDEFAULTFOLD;
}

int verify_sicssq_reproducibility(int fold, int N, float complex* X, int incX, float refscl, float refssq, float_indexed *iref, int max_num_blocks) {
  // GENERATE DATA
  int i;
  float resscl = 0.0;
  float resssq;
  float_indexed *ires = idxd_sialloc(fold);
  int num_blocks = 1;

  int block_N = (N + num_blocks - 1) / num_blocks;

  num_blocks = 1;
  while (num_blocks < N && num_blocks <= max_num_blocks) {
    idxd_sisetzero(fold, ires);
    if (num_blocks == 1)
      resscl = idxdBLAS_sicssq(fold, N, X, incX, resscl, ires);
    else {
      block_N = (N + num_blocks - 1) / num_blocks;
      for (i = 0; i < N; i += block_N) {
        block_N = block_N < N - i ? block_N : (N-i);
        resscl = idxdBLAS_sicssq(fold, block_N, X + i * incX, incX, resscl, ires);
      }
    }
    resssq = idxd_ssiconv(fold, ires);
    if (resssq != refssq || resscl != refscl) {
      printf("idxdBLAS_sicssq(X)[num_blocks=%d,block_N=%d] = %g * (%g) != %g * (%g)\n", num_blocks, block_N, resscl, resssq, refscl, refssq);
      printf("ref float_indexed:\n");
      idxd_siprint(fold, iref);
      printf("\nres float_indexed:\n");
      idxd_siprint(fold, ires);
      printf("\n");
      return 1;
    }
    num_blocks *= 2;
  }
  free(ires);
  return 0;
}

int vecvec_fill_show_help(void){
  verify_sicssq_options_initialize();

  opt_show_option(max_blocks);
  opt_show_option(shuffles);
  opt_show_option(fold);
  return 0;
}

const char* vecvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  verify_sicssq_options_initialize();

  opt_eval_option(argc, argv, &fold);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Verify sicssq reproducibility fold=%d", fold._int.value);
  return name_buffer;
}

int vecvec_fill_test(int argc, char** argv, int N, int FillX, double RealScaleX, double ImagScaleX, int incX, int FillY, double RealScaleY, double ImagScaleY, int incY){
  int rc = 0;
  int i;
  float refscl = 0.0;
  float refssq;
  float_indexed *iref;
  int max_num_blocks;

  verify_sicssq_options_initialize();

  util_random_seed();

  opt_eval_option(argc, argv, &max_blocks);
  opt_eval_option(argc, argv, &shuffles);
  opt_eval_option(argc, argv, &fold);

  iref = idxd_sialloc(fold._int.value);

  float complex *X = util_cvec_alloc(N, incX);
  float complex *Y = util_cvec_alloc(N, incY);
  int *P;


  max_num_blocks = max_blocks._int.value;

  util_cvec_fill(N, X, incX, FillX, RealScaleX, ImagScaleX);
  util_cvec_fill(N, Y, incY, FillY, RealScaleY, ImagScaleY);


  //compute with unpermuted data
  idxd_sisetzero(fold._int.value, iref);
  refscl = idxdBLAS_sicssq(fold._int.value, N, X, incX, refscl, iref);
  refssq = idxd_ssiconv(fold._int.value, iref);

  P = util_identity_permutation(N);
  util_cvec_reverse(N, X, incX, P, 1);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_sicssq_reproducibility(fold._int.value, N, X, incX, refscl, refssq, iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_sort(N, X, incX, P, 1, util_Increasing);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_sicssq_reproducibility(fold._int.value, N, X, incX, refscl, refssq, iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_sort(N, X, incX, P, 1, util_Decreasing);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_sicssq_reproducibility(fold._int.value, N, X, incX, refscl, refssq, iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_sort(N, X, incX, P, 1, util_Increasing_Magnitude);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_sicssq_reproducibility(fold._int.value, N, X, incX, refscl, refssq, iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(N);
  util_cvec_sort(N, X, incX, P, 1, util_Decreasing_Magnitude);
  util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
  free(P);

  rc = verify_sicssq_reproducibility(fold._int.value, N, X, incX, refscl, refssq, iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  for(i = 0; i < shuffles._int.value; i++){
    P = util_identity_permutation(N);
    util_cvec_shuffle(N, X, incX, P, 1);
    util_cvec_permute(N, Y, incY, P, 1, NULL, 1);
    free(P);

    rc = verify_sicssq_reproducibility(fold._int.value, N, X, incX, refscl, refssq, iref, max_num_blocks);
    if(rc != 0){
      return rc;
    }
  }

  free(iref);
  free(X);
  free(Y);

  return rc;
}
