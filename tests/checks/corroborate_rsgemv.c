#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <binned.h>
#include <binnedBLAS.h>
#include <reproBLAS.h>

#include "../common/test_opt.h"
#include "../common/test_matvec_fill_header.h"

#include "wrap_rsgemv.h"

static opt_option max_blocks;
static opt_option shuffles;
static opt_option fold;

static void corroborate_rsgemv_options_initialize(void){
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
  fold._int.max               = binned_SBMAXFOLD;
  fold._int.value             = SIDEFAULTFOLD;
}

int corroborate_rsgemv(int fold, char Order, char TransA, int M, int N, float alpha, float *A, int lda, float* X, int incX, float beta, float *Y, float_binned *YI, int incY, float *ref, int max_num_blocks) {

  int i;
  int num_blocks = 1;
  int block_opN;
  int opN;
  int opM;
  switch(TransA){
    case 'n':
    case 'N':
      opN = N;
      opM = M;
    break;
    default:
      opN = M;
      opM = N;
    break;
  }

  float *res = util_svec_alloc(opM, incY);
  float_binned *Ires = malloc(opM * incY * binned_sbsbze(fold));

  num_blocks = 1;
  while (num_blocks < opN && num_blocks <= max_num_blocks) {
    memcpy(res, Y, opM * incY * sizeof(float));
    memcpy(Ires, YI, opM * incY * binned_sbsbze(fold));
    if (num_blocks == 1){
      wrap_rsgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, beta, res, incY);
    }else {
      switch(TransA){
        case 'n':
        case 'N':
          {
            block_opN = (N + num_blocks - 1) / num_blocks;
            for (i = 0; i < N; i += block_opN) {
              block_opN = block_opN < N - i ? block_opN : (N-i);
              switch(Order){
                case 'r':
                case 'R':
                  binnedBLAS_sbsgemv(fold, Order, TransA, M, block_opN, alpha, A + i, lda, X + i * incX, incX, Ires, incY);
                  break;
                default:
                  binnedBLAS_sbsgemv(fold, Order, TransA, M, block_opN, alpha, A + i * lda, lda, X + i * incX, incX, Ires, incY);
                  break;
              }
            }
          }
          break;
        default:
          {
            block_opN = (M + num_blocks - 1) / num_blocks;
            for (i = 0; i < M; i += block_opN) {
              block_opN = block_opN < M - i ? block_opN : (M-i);
              switch(Order){
                case 'r':
                case 'R':
                  binnedBLAS_sbsgemv(fold, Order, TransA, block_opN, N, alpha, A + i * lda, lda, X + i * incX, incX, Ires, incY);
                  break;
                default:
                  binnedBLAS_sbsgemv(fold, Order, TransA, block_opN, N, alpha, A + i, lda, X + i * incX, incX, Ires, incY);
                  break;
              }
            }
          }
          break;
      }
      for(i = 0; i < opM; i++){
        res[i * incY] = binned_ssbconv(fold, Ires + i * incY * binned_sbnum(fold));
      }
    }
    for(i = 0; i < opM; i++){
      if(res[i * incY] != ref[i * incY]){
        printf("reproBLAS_rsgemv(A, X, Y)[num_blocks=%d] = %g != %g\n", num_blocks, res[i * incY], ref[i * incY]);
        return 1;
      }
    }
    num_blocks *= 2;
  }
  free(res);
  free(Ires);
  return 0;
}

int matvec_fill_show_help(void){
  corroborate_rsgemv_options_initialize();

  opt_show_option(fold);
  opt_show_option(max_blocks);
  opt_show_option(shuffles);
  return 0;
}

const char* matvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  corroborate_rsgemv_options_initialize();

  opt_eval_option(argc, argv, &fold);
  opt_eval_option(argc, argv, &max_blocks);

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Corroborate rsgemv results fold=%d", fold._int.value);
  return name_buffer;
}

int matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillX, double RealScaleX, double ImagScaleX, int incX, double RealBeta, double ImagBeta, int FillY, double RealScaleY, double ImagScaleY, int incY){
  (void)ImagAlpha;
  (void)ImagBeta;
  int rc = 0;
  int i;

  corroborate_rsgemv_options_initialize();

  opt_eval_option(argc, argv, &fold);
  opt_eval_option(argc, argv, &max_blocks);
  opt_eval_option(argc, argv, &shuffles);

  util_random_seed();
  int opN;
  int opM;
  char NTransA;
  switch(TransA){
    case 'n':
    case 'N':
      opN = N;
      opM = M;
      NTransA = 'T';
    break;
    default:
      opN = M;
      opM = N;
      NTransA = 'N';
    break;
  }

  float *A  = util_smat_alloc(Order, M, N, lda);
  float *X  = util_svec_alloc(opN, incX);
  float *Y  = util_svec_alloc(opM, incY);
  float_binned *YI = (float_binned*)malloc(opM * incY * binned_sbsbze(fold._int.value));

  int *P;

  util_smat_fill(Order, 'n', M, N, A, lda, FillA, RealScaleA, ImagScaleA);
  util_svec_fill(opN, X, incX, FillX, RealScaleX, ImagScaleX);
  util_svec_fill(opM, Y, incY, FillY, RealScaleY, ImagScaleY);
  if(RealBeta == 0.0){
    memset(YI, 0, opM * binned_sbsbze(fold._int.value));
  }else{
    for(i = 0; i < opM; i++){
      binned_sbsconv(fold._int.value, Y[i * incY] * RealBeta, YI + i * incY * binned_sbnum(fold._int.value));
    }
  }
  float *ref  = (float*)malloc(opM * incY * sizeof(float));

  //compute with unpermuted data
  memcpy(ref, Y, opM * incY * sizeof(float));

  wrap_ref_rsgemv(fold._int.value, Order, TransA, M, N, RealAlpha, A, lda, X, incX, RealBeta, ref, incY);

  rc = corroborate_rsgemv(fold._int.value, Order, TransA, M, N, RealAlpha, A, lda, X, incX, RealBeta, Y, YI, incY, ref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(opN);
  util_svec_reverse(opN, X, incX, P, 1);
  util_smat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = corroborate_rsgemv(fold._int.value, Order, TransA, M, N, RealAlpha, A, lda, X, incX, RealBeta, Y, YI, incY, ref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(opN);
  util_svec_sort(opN, X, incX, P, 1, util_Increasing);
  util_smat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = corroborate_rsgemv(fold._int.value, Order, TransA, M, N, RealAlpha, A, lda, X, incX, RealBeta, Y, YI, incY, ref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(opN);
  util_svec_sort(opN, X, incX, P, 1, util_Decreasing);
  util_smat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = corroborate_rsgemv(fold._int.value, Order, TransA, M, N, RealAlpha, A, lda, X, incX, RealBeta, Y, YI, incY, ref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(opN);
  util_svec_sort(opN, X, incX, P, 1, util_Increasing_Magnitude);
  util_smat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = corroborate_rsgemv(fold._int.value, Order, TransA, M, N, RealAlpha, A, lda, X, incX, RealBeta, Y, YI, incY, ref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(opN);
  util_svec_sort(opN, X, incX, P, 1, util_Decreasing_Magnitude);
  util_smat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = corroborate_rsgemv(fold._int.value, Order, TransA, M, N, RealAlpha, A, lda, X, incX, RealBeta, Y, YI, incY, ref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  for(i = 0; i < shuffles._int.value; i++){
    P = util_identity_permutation(opN);
    util_svec_shuffle(opN, X, incX, P, 1);
    util_smat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
    free(P);

    rc = corroborate_rsgemv(fold._int.value, Order, TransA, M, N, RealAlpha, A, lda, X, incX, RealBeta, Y, YI, incY, ref, max_blocks._int.value);
    if(rc != 0){
      return rc;
    }
  }

  free(A);
  free(X);
  free(Y);
  free(YI);
  free(ref);

  return rc;
}
