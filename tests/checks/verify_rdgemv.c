#include <indexedBLAS.h>
#include <indexed.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../common/test_opt.h"

#include "../common/test_matvec_fill_header.h"

static opt_option incY;
static opt_option FillY;
static opt_option ScaleY;
static opt_option CondY;
static opt_option alpha;
static opt_option beta;

static void verify_rdgemv_options_initialize(void){
  incY._int.header.type       = opt_int;
  incY._int.header.short_name = 'y';
  incY._int.header.long_name  = "incY";
  incY._int.header.help       = "Y vector increment";
  incY._int.required          = 0;
  incY._int.min               = 1;
  incY._int.max               = INT_MAX;
  incY._int.value             = 1;

  FillY._named.header.type       = opt_named;
  FillY._named.header.short_name = 'j';
  FillY._named.header.long_name  = "FillY";
  FillY._named.header.help       = "Y fill type";
  FillY._named.required          = 0;
  FillY._named.n_names           = (int)util_vec_fill_n_names;
  FillY._named.names             = (char**)util_vec_fill_names;
  FillY._named.descs             = (char**)util_vec_fill_descs;
  FillY._named.value             = 0;

  ScaleY._double.header.type       = opt_double;
  ScaleY._double.header.short_name = 'v';
  ScaleY._double.header.long_name  = "ScaleY";
  ScaleY._double.header.help       = "Y scale";
  ScaleY._double.required          = 0;
  ScaleY._double.min               = 0;
  ScaleY._double.max               = DBL_MAX;
  ScaleY._double.value             = 1.0;

  CondY._double.header.type       = opt_double;
  CondY._double.header.short_name = 'e';
  CondY._double.header.long_name  = "CondY";
  CondY._double.header.help       = "Y condition number";
  CondY._double.required          = 0;
  CondY._double.min               = 1.0;
  CondY._double.max               = DBL_MAX;
  CondY._double.value             = 1e3;

  alpha._double.header.type       = opt_double;
  alpha._double.header.short_name = 'l';
  alpha._double.header.long_name  = "alpha";
  alpha._double.header.help       = "alpha";
  alpha._double.required          = 0;
  alpha._double.min               = 1.0;
  alpha._double.max               = DBL_MAX;
  alpha._double.value             = 1e3;

  beta._double.header.type       = opt_double;
  beta._double.header.short_name = 'm';
  beta._double.header.long_name  = "beta";
  beta._double.header.help       = "beta";
  beta._double.required          = 0;
  beta._double.min               = 1.0;
  beta._double.max               = DBL_MAX;
  beta._double.value             = 1e3;
}

void wrap_rdgemv(const char Order,
                 const char TransA, const int M, const int N,
                 const double *A, const double alpha, const int lda,
                 const double *X, const int incX,
                 const double beta, double *Y, const int incY){
  (void)alpha;
  (void)beta;
  rblas_order_t o;
  rblas_transpose_t t;
  switch(Order){
    case 'R':
      o = rblas_Row_Major;
      break;
    default:
      o = rblas_Col_Major;
      break;
  }
  switch(TransA){
    case 'N':
      t = rblas_No_Trans;
      break;
    default:
      t = rblas_Trans;
      break;
  }
  rdgemv(o, t, M, N, A, lda, X, incX, Y, incY);
}

void wrap_dgemvI(const char Order,
                 const char TransA, const int M, const int N,
                 const double *A, const double alpha, const int lda,
                 const double *X, const int incX,
                 const double beta, double_indexed *Y, const int incY){
  (void)alpha;
  (void)beta;
  rblas_order_t o;
  rblas_transpose_t t;
  switch(Order){
    case 'R':
      o = rblas_Row_Major;
      break;
    default:
      o = rblas_Col_Major;
      break;
  }
  switch(TransA){
    case 'N':
      t = rblas_No_Trans;
      break;
    default:
      t = rblas_Trans;
      break;
  }
  dgemvI(DEFAULT_FOLD, o, t, M, N, A, lda, X, incX, Y, incY);
}

int verify_dgemv_reproducibility(char Order, char TransA, int M, int N, int NX, int NY, double alpha, double *A, int lda, double* X, int incX, double beta, double *Y, double_indexed *YI, int incY, double *ref, double_indexed *Iref, int max_num_blocks) {
  (void)NX;
  (void)alpha;
  (void)beta;

  // GENERATE DATA
  int i;
  int num_blocks = 1;
  int block_N;

  double *res = malloc(NY * incY * sizeof(double));
  double_indexed *Ires = malloc(NY * incY * disize(DEFAULT_FOLD));

  num_blocks = 1;
  while (num_blocks < N && num_blocks <= max_num_blocks) {
    //compute with unpermuted data
    memcpy(res, Y, NY * incY * sizeof(double));
    memcpy(Ires, YI, NY * incY * disize(DEFAULT_FOLD));
    if (num_blocks == 1)
      wrap_rdgemv(Order, TransA, M, N, A, alpha, lda, X, incX, beta, res, incY);
    else {
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
                  wrap_dgemvI(Order, TransA, M, block_N, A + i, alpha, lda, X + i * incX, incX, beta, Ires, incY);
                  break;
                default:
                  wrap_dgemvI(Order, TransA, M, block_N, A + lda * i, alpha, lda, X + i * incX, incX, beta, Ires, incY);
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
                  wrap_dgemvI(Order, TransA, block_N, N, A + i * lda, alpha, lda, X + i * incX, incX, beta, Ires, incY);
                  break;
                default:
                  wrap_dgemvI(Order, TransA, block_N, N, A + i, alpha, lda, X + i * incX, incX, beta, Ires, incY);
                  break;
              }
            }
          }
          break;
      }
      for(i = 0; i < NY; i++){
        res[i * incY] = ddiconv(DEFAULT_FOLD, Ires + i * incY * dinum(DEFAULT_FOLD));
      }
      for(i = 0; i < NY; i++){
        if(res[i * incY] != ref[i * incY]){
          printf("dgemv(A, X, Y)[num_blocks=%d,block_N=%d] = %g != %g\n", num_blocks, block_N, res[i * incY], ref[i * incY]);
          if (num_blocks != 1) {
            printf("Ref I_double:\n");
            diprint(DEFAULT_FOLD, Iref + i * incY * dinum(DEFAULT_FOLD));
            printf("\nRes I_double:\n");
            diprint(DEFAULT_FOLD, Ires + i * incY * dinum(DEFAULT_FOLD));
            printf("\n");
          }
          return 1;
        }
      }
    }
    num_blocks *= 2;
  }
  return 0;
}

int matvec_fill_show_help(void){
  verify_rdgemv_options_initialize();

  opt_show_option(incY);
  opt_show_option(FillY);
  opt_show_option(ScaleY);
  opt_show_option(CondY);
  opt_show_option(alpha);
  opt_show_option(beta);
  return 0;
}

const char* matvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  verify_rdgemv_options_initialize();

  opt_eval_option(argc, argv, &incY);
  opt_eval_option(argc, argv, &FillY);
  opt_eval_option(argc, argv, &ScaleY);
  opt_eval_option(argc, argv, &CondY);
  opt_eval_option(argc, argv, &alpha);
  opt_eval_option(argc, argv, &beta);

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Verify dgemv reproducibility");
  return name_buffer;
}

int matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, int FillA, double ScaleA, double CondA, int lda, int FillX, double ScaleX, double CondX, int incX){
  int rc = 0;
  int i;
  int max_num_blocks = 1024;

  verify_rdgemv_options_initialize();

  opt_eval_option(argc, argv, &incY);
  opt_eval_option(argc, argv, &FillY);
  opt_eval_option(argc, argv, &ScaleY);
  opt_eval_option(argc, argv, &CondY);
  opt_eval_option(argc, argv, &alpha);
  opt_eval_option(argc, argv, &beta);

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

  double *A  = util_dmat_alloc(Order, M, N, lda);
  double *X  = util_dvec_alloc(NX, incX);
  double *Y  = util_dvec_alloc(NY, incY._int.value);
  double_indexed *YI = (double_indexed*)malloc(NY * incY._int.value * disize(DEFAULT_FOLD));

  int *P;

  util_dmat_fill(Order, 'n', M, N, A, lda, FillA, ScaleA, CondA);
  util_dvec_fill(NX, X, incX, FillX, ScaleX, CondX);
  util_dvec_fill(NY, Y, incY._int.value, FillY._named.value, ScaleY._double.value, CondY._double.value);
  for(i = 0; i < NY; i++){
    didconv(DEFAULT_FOLD, Y[i * incY._int.value], YI + i * incY._int.value * dinum(DEFAULT_FOLD));
  }
  double *ref  = (double*)malloc(NY * incY._int.value * sizeof(double));
  double_indexed *Iref = (double_indexed*)malloc(NY * incY._int.value * disize(DEFAULT_FOLD));

  //compute with unpermuted data
  memcpy(ref, Y, NY * incY._int.value * sizeof(double));
  memcpy(Iref, YI, NY * incY._int.value * disize(DEFAULT_FOLD));

  wrap_rdgemv(Order, TransA, M, N, A, alpha._double.value, lda, X, incX, beta._double.value, ref, incY._int.value);
  wrap_dgemvI(Order, TransA, M, N, A, alpha._double.value, lda, X, incX, beta._double.value, Iref, incY._int.value);

  P = util_identity_permutation(NX);
  util_dvec_reverse(NX, X, incX, P, 1);
  util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_dgemv_reproducibility(Order, TransA, M, N, NX, NY, alpha._double.value, A, lda, X, incX, beta._double.value, Y, YI, incY._int.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_dvec_sort(NX, X, incX, P, 1, util_Increasing);
  util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_dgemv_reproducibility(Order, TransA, M, N, NX, NY, alpha._double.value, A, lda, X, incX, beta._double.value, Y, YI, incY._int.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_dvec_sort(NX, X, incX, P, 1, util_Decreasing);
  util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_dgemv_reproducibility(Order, TransA, M, N, NX, NY, alpha._double.value, A, lda, X, incX, beta._double.value, Y, YI, incY._int.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_dvec_sort(NX, X, incX, P, 1, util_Increasing_Magnitude);
  util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_dgemv_reproducibility(Order, TransA, M, N, NX, NY, alpha._double.value, A, lda, X, incX, beta._double.value, Y, YI, incY._int.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_dvec_sort(NX, X, incX, P, 1, util_Decreasing_Magnitude);
  util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_dgemv_reproducibility(Order, TransA, M, N, NX, NY, alpha._double.value, A, lda, X, incX, beta._double.value, Y, YI, incY._int.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_dvec_shuffle(NX, X, incX, P, 1);
  util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_dgemv_reproducibility(Order, TransA, M, N, NX, NY, alpha._double.value, A, lda, X, incX, beta._double.value, Y, YI, incY._int.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_dvec_shuffle(NX, X, incX, P, 1);
  util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_dgemv_reproducibility(Order, TransA, M, N, NX, NY, alpha._double.value, A, lda, X, incX, beta._double.value, Y, YI, incY._int.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_dvec_shuffle(NX, X, incX, P, 1);
  util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_dgemv_reproducibility(Order, TransA, M, N, NX, NY, alpha._double.value, A, lda, X, incX, beta._double.value, Y, YI, incY._int.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  P = util_identity_permutation(NX);
  util_dvec_shuffle(NX, X, incX, P, 1);
  util_dmat_row_permute(Order, NTransA, M, N, A, lda, P, 1, NULL, 1);
  free(P);

  rc = verify_dgemv_reproducibility(Order, TransA, M, N, NX, NY, alpha._double.value, A, lda, X, incX, beta._double.value, Y, YI, incY._int.value, ref, Iref, max_num_blocks);
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
