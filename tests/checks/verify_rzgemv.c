#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <indexed.h>
#include <indexedBLAS.h>
#include <reproBLAS.h>

#include "../common/test_opt.h"
#include "../common/test_matvec_fill_header.h"

#include "wrap_rzgemv.h"

static opt_option max_blocks;
static opt_option fold;

static void verify_rzgemv_options_initialize(void){
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
  fold._int.max               = DIMAXFOLD;
  fold._int.value             = DIDEFAULTFOLD;
}

int verify_zgemv_reproducibility(int fold, char Order, char TransA, int M, int N, int NX, int NY, double complex *alpha, double complex *A, int lda, double complex *X, int incX, double complex *beta, double complex *Y, double_complex_indexed *YI, int incY, double complex *ref, double_complex_indexed *Iref, int max_num_blocks) {
  int i;
  int num_blocks = 1;
  int block_N;

  double complex *res = malloc(NY * incY * sizeof(double complex));
  double_complex_indexed *Ires = malloc(NY * incY * zisize(fold));

  num_blocks = 1;
  while (num_blocks < N && num_blocks <= max_num_blocks) {
    memcpy(res, Y, NY * incY * sizeof(double complex));
    memcpy(Ires, YI, NY * incY * zisize(fold));
    if (num_blocks == 1){
      wrap_rzgemv(fold, Order, TransA, M, N, alpha, A, lda, X, incX, beta, res, incY);
    }else {
      switch(TransA){
        case 'n':
        case 'N':
          {
            block_N = (N + num_blocks - 1) / num_blocks;
            for (i = 0; i < N; i += block_N) {
              block_N = block_N < N - i ? block_N : (N-i);
              switch(Order){
                case 'r':
                case 'R':
                  zizgemv(fold, Order, TransA, M, block_N, (void*)alpha, (void*)(A + i), lda, (void*)(X + i * incX), incX, Ires, incY);
                  printf("i %d p1 %g\n", i, Ires[3]);
                  break;
                default:
                  zizgemv(fold, Order, TransA, M, block_N, (void*)alpha, (void*)(A + i * lda), lda, (void*)(X + i * incX), incX, Ires, incY);
                  break;
              }
            }
          }
          break;
        default:
          {
            block_N = (M + num_blocks - 1) / num_blocks;
            for (i = 0; i < M; i += block_N) {
              block_N = block_N < M - i ? block_N : (M-i);
              switch(Order){
                case 'r':
                case 'R':
                  zizgemv(fold, Order, TransA, block_N, N, (void*)alpha, (void*)(A + i * lda), lda, (void*)(X + i * incX), incX, Ires, incY);
                  break;
                default:
                  zizgemv(fold, Order, TransA, block_N, N, (void*)alpha, (void*)(A + i), lda, (void*)(X + i * incX), incX, Ires, incY);
                  break;
              }
            }
          }
          break;
      }
      for(i = 0; i < NY; i++){
        zziconv_sub(fold, Ires + i * incY * zinum(fold), (void*)(res + i * incY));
      }
    }
    for(i = 0; i < NY; i++){
      if(res[i * incY] != ref[i * incY]){
        printf("rzgemv(A, X, Y)[num_blocks=%d,block_N=%d] = %g + %gi != %g + %gi\n", num_blocks, block_N, creal(res[i * incY]), cimag(res[i * incY]), creal(ref[i * incY]), cimag(ref[i * incY]));
        if (num_blocks != 1) {
          printf("Ref I_double_complex:\n");
          ziprint(fold, Iref + i * incY * zinum(fold));
          printf("\nRes I_double_complex:\n");
          ziprint(fold, Ires + i * incY * zinum(fold));
          printf("\n");
        }
        return 1;
      }
    }
    num_blocks *= 2;
  }
  return 0;
}

int matvec_fill_show_help(void){
  verify_rzgemv_options_initialize();

  opt_show_option(fold);
  opt_show_option(max_blocks);
  return 0;
}

const char* matvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  verify_rzgemv_options_initialize();

  opt_eval_option(argc, argv, &fold);
  opt_eval_option(argc, argv, &max_blocks);

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Verify rzgemv reproducibility fold=%d", fold._int.value);
  return name_buffer;
}

int matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillX, double RealScaleX, double ImagScaleX, int incX, double RealBeta, double ImagBeta, int FillY, double RealScaleY, double ImagScaleY, int incY){
  int rc = 0;
  int i;

  verify_rzgemv_options_initialize();

  opt_eval_option(argc, argv, &fold);
  opt_eval_option(argc, argv, &max_blocks);

  util_random_seed();
  int NX;
  int NY;
  char NTransA;
  switch(TransA){
    case 'n':
    case 'N':
      NX = N;
      NY = M;
      NTransA = 'T';
    break;
    default:
      NX = M;
      NY = N;
      NTransA = 'N';
    break;
  }

  double complex *A = util_zmat_alloc(Order, M, N, lda);
  double complex *X = util_zvec_alloc(NX, incX);
  double complex *Y = util_zvec_alloc(NY, incY);
  double complex alpha = RealAlpha + I * ImagAlpha;
  double complex beta = RealBeta + I * ImagBeta;
  double complex betaY;
  double_complex_indexed *YI = (double_complex_indexed*)malloc(NY * incY * zisize(fold._int.value));

  int *P;

  util_zmat_fill(Order, 'n', M, N, A, lda, FillA, RealScaleA, ImagScaleA);
  util_zvec_fill(NX, X, incX, FillX, RealScaleX, ImagScaleX);
  util_zvec_fill(NY, Y, incY, FillY, RealScaleY, ImagScaleY);
  for(i = 0; i < NY; i++){
    betaY = Y[i * incY] * beta;
    zizconv(fold._int.value, (void*)&betaY, YI + i * incY * zinum(fold._int.value));
  }
  double complex *ref  = (double complex*)malloc(NY * incY * sizeof(double complex));
  double_complex_indexed *Iref = (double_complex_indexed*)malloc(NY * incY * zisize(fold._int.value));

  //compute with unpermuted data
  memcpy(ref, Y, NY * incY * sizeof(double complex));
  memcpy(Iref, YI, NY * incY * zisize(fold._int.value));

  wrap_rzgemv(fold._int.value, Order, TransA, M, N, &alpha, A, lda, X, incX, &beta, ref, incY);
  zizgemv(fold._int.value, Order, TransA, M, N, &alpha, A, lda, X, incX, Iref, incY);

  P = util_identity_permutation(NX);
  util_zvec_reverse(NX, X, incX, P, 1);
  util_zmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_zgemv_reproducibility(fold._int.value, Order, TransA, M, N, NX, NY, &alpha, A, lda, X, incX, &beta, Y, YI, incY, ref, Iref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_zvec_sort(NX, X, incX, P, 1, util_Increasing);
  util_zmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_zgemv_reproducibility(fold._int.value, Order, TransA, M, N, NX, NY, &alpha, A, lda, X, incX, &beta, Y, YI, incY, ref, Iref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_zvec_sort(NX, X, incX, P, 1, util_Decreasing);
  util_zmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_zgemv_reproducibility(fold._int.value, Order, TransA, M, N, NX, NY, &alpha, A, lda, X, incX, &beta, Y, YI, incY, ref, Iref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_zvec_sort(NX, X, incX, P, 1, util_Increasing_Magnitude);
  util_zmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_zgemv_reproducibility(fold._int.value, Order, TransA, M, N, NX, NY, &alpha, A, lda, X, incX, &beta, Y, YI, incY, ref, Iref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_zvec_sort(NX, X, incX, P, 1, util_Decreasing_Magnitude);
  util_zmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_zgemv_reproducibility(fold._int.value, Order, TransA, M, N, NX, NY, &alpha, A, lda, X, incX, &beta, Y, YI, incY, ref, Iref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_zvec_shuffle(NX, X, incX, P, 1);
  util_zmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_zgemv_reproducibility(fold._int.value, Order, TransA, M, N, NX, NY, &alpha, A, lda, X, incX, &beta, Y, YI, incY, ref, Iref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_zvec_shuffle(NX, X, incX, P, 1);
  util_zmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_zgemv_reproducibility(fold._int.value, Order, TransA, M, N, NX, NY, &alpha, A, lda, X, incX, &beta, Y, YI, incY, ref, Iref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_zvec_shuffle(NX, X, incX, P, 1);
  util_zmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_zgemv_reproducibility(fold._int.value, Order, TransA, M, N, NX, NY, &alpha, A, lda, X, incX, &beta, Y, YI, incY, ref, Iref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_zvec_shuffle(NX, X, incX, P, 1);
  util_zmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_zgemv_reproducibility(fold._int.value, Order, TransA, M, N, NX, NY, &alpha, A, lda, X, incX, &beta, Y, YI, incY, ref, Iref, max_blocks._int.value);
  if(rc != 0){
    return rc;
  }

  free(A);
  free(X);
  free(Y);
  free(ref);
  free(Iref);

  return rc;
}
