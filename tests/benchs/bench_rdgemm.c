#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <idxdBLAS.h>
#include <reproBLAS.h>

#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_metric.h"

#include "../../config.h"

#include "bench_matmat_fill_header.h"

static opt_option fold;

static void bench_rdgemm_options_initialize(void){
  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 0;
  fold._int.max               = DIMAXFOLD;
  fold._int.value             = DIDEFAULTFOLD;
}

int bench_matmat_fill_show_help(void){
  bench_rdgemm_options_initialize();

  opt_show_option(fold);
  return 0;
}

const char* bench_matmat_fill_name(int argc, char** argv){
  (void)argc;
  (void)argv;
  bench_rdgemm_options_initialize();

  opt_eval_option(argc, argv, &fold);

  static char name_buffer[MAX_LINE];
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [rdgemm] (fold = %d)", fold._int.value);
  return name_buffer;
}

int bench_matmat_fill_test(int argc, char** argv, char Order, char TransA, char TransB, int M, int N, int K, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillB, double RealScaleB, double ImagScaleB, int ldb, double RealBeta, double ImagBeta, int FillC, double RealScaleC, double ImagScaleC, int ldc, int trials){
  int rc = 0;
  int t;
  int i;
  int j;
  int k;

  bench_rdgemm_options_initialize();

  opt_eval_option(argc, argv, &fold);

  util_random_seed();

  double *A  = util_dmat_alloc(Order, M, K, lda);
  double *B  = util_dmat_alloc(Order, K, N, ldb);
  double *C  = util_dmat_alloc(Order, M, N, ldc);
  double *res  = util_dmat_alloc(Order, M, N, ldc);
  double alpha = RealAlpha;
  double beta = RealBeta;
  double_indexed *CI;

  util_dmat_fill(Order, TransA, M, N, A, lda, FillA, RealScaleA, ImagScaleA);
  util_dmat_fill(Order, TransB, M, N, A, ldb, FillB, RealScaleB, ImagScaleB);
  util_dmat_fill(Order, 'n', M, N, A, ldc, FillC, RealScaleC, ImagScaleC);

  if(fold._int.value == DIDEFAULTFOLD){
    for(i = 0; i < trials; i++){
      switch(Order){
        case 'r':
        case 'R':
          memcpy(res, C, M * ldc * sizeof(double));
          break;
        default:
          memcpy(res, C, ldc * N * sizeof(double));
          break;
      }
      time_tic();
      reproBLAS_rdgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, res, ldc);
      time_toc();
    }
  }else{
    for(t = 0; t < trials; t++){
      switch(Order){
        case 'r':
        case 'R':
          memcpy(res, C, M * ldc * sizeof(double));
          break;
        default:
          memcpy(res, C, ldc * N * sizeof(double));
          break;
      }
      time_tic();
      CI = (double_indexed*)malloc(M * N * idxd_disize(DIDEFAULTFOLD));
      switch(Order){
        case 'r':
        case 'R':
          if(beta == 1.0){
            for(i = 0; i < M; i++){
              for(j = 0; j < N; j++){
                idxd_didconv(fold._int.value, res[i * ldc + j], CI + (i * N + j) * idxd_dinum(fold._int.value));
              }
            }
          }else if(beta == 0.0){
            memset(CI, 0.0, M * N * idxd_disize(fold._int.value));
          }else{
            for(i = 0; i < M; i++){
              for(j = 0; j < N; j++){
                idxd_didconv(fold._int.value, res[i * ldc + j] * beta, CI + (i * N + j) * idxd_dinum(fold._int.value));
              }
            }
          }
          idxdBLAS_didgemm(fold._int.value, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, CI, N);
          for(i = 0; i < M; i++){
            for(j = 0; j < N; j++){
              res[i * ldc + j] = idxd_ddiconv(fold._int.value, CI + (i * N + j) * idxd_dinum(fold._int.value));
            }
          }
          break;
        default:
          if(beta == 1.0){
            for(j = 0; j < N; j++){
              for(i = 0; i < M; i++){
                idxd_didconv(fold._int.value, res[j * ldc + i], CI + (j * M + i) * idxd_dinum(fold._int.value));
              }
            }
          }else if(beta == 0.0){
            memset(CI, 0.0, M * N * idxd_disize(fold._int.value));
          }else{
            for(j = 0; j < N; j++){
              for(i = 0; i < M; i++){
                idxd_didconv(fold._int.value, res[j * ldc + i] * beta, CI + (j * M + i) * idxd_dinum(fold._int.value));
              }
            }
          }
          idxdBLAS_didgemm(fold._int.value, Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, CI, M);
          for(j = 0; j < N; j++){
            for(i = 0; i < M; i++){
              res[j * ldc + i] = idxd_ddiconv(fold._int.value, CI + (j * M + i) * idxd_dinum(fold._int.value));
            }
          }
          break;
      }
      free(CI);
      time_toc();
    }
  }

  double dM = (double)M;
  double dN = (double)N;
  double dK = (double)K;
  metric_load_double("time", time_read());
  metric_load_double("trials", (double)(trials));
  metric_load_double("input", (double)(dM * dK + dK * dN + dM * dN));
  metric_load_double("output", (double)dN * dM);
  metric_load_double("d_mul", (double)(dN * dM * dK));
  metric_load_double("d_add", (double)((3 * fold._int.value - 2) * dN * dM * dK));
  metric_load_double("d_orb", (double)(fold._int.value * dN * dM * dK));
  metric_dump();

  free(A);
  free(B);
  free(C);
  free(res);
  return rc;
}
