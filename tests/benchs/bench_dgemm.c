#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <idxdBLAS.h>
#include <reproBLAS.h>

#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_metric.h"
#include "../common/test_BLAS.h"

#include "../../config.h"

#include "bench_matmat_fill_header.h"

int bench_matmat_fill_show_help(void){
  return 0;
}

const char* bench_matmat_fill_name(int argc, char** argv){
  (void)argc;
  (void)argv;

  static char name_buffer[MAX_LINE];
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [rdgemm]");
  return name_buffer;
}

int bench_matmat_fill_test(int argc, char** argv, char Order, char TransA, char TransB, int M, int N, int K, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillB, double RealScaleB, double ImagScaleB, int ldb, double RealBeta, double ImagBeta, int FillC, double RealScaleC, double ImagScaleC, int ldc, int trials){
  int rc = 0;
  int i;

  util_random_seed();

  double *A  = util_dmat_alloc(Order, M, K, lda);
  double *B  = util_dmat_alloc(Order, K, N, ldb);
  double *C  = util_dmat_alloc(Order, M, N, ldc);
  double *res  = util_dmat_alloc(Order, M, N, ldc);
  double alpha = RealAlpha;
  double beta = RealBeta;

  util_dmat_fill(Order, TransA, M, N, A, lda, FillA, RealScaleA, ImagScaleA);
  util_dmat_fill(Order, TransB, M, N, B, ldb, FillB, RealScaleB, ImagScaleB);
  util_dmat_fill(Order, 'n', M, N, C, ldc, FillC, RealScaleC, ImagScaleC);

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
    CALL_DGEMM(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, res, ldc);
    time_toc();
  }

  double dM = (double)M;
  double dN = (double)N;
  double dK = (double)K;
  metric_load_double("time", time_read());
  metric_load_double("trials", (double)(trials));
  metric_load_double("input", (double)(dM * dK + dK * dN + dM * dN));
  metric_load_double("output", (double)dN * dM);
  metric_load_double("d_fma", (double)(dN * dM * dK));
  metric_dump();

  free(A);
  free(B);
  free(C);
  free(res);
  return rc;
}
