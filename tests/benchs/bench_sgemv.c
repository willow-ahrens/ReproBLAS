#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_metric.h"
#include "../common/test_BLAS.h"

#include "../../config.h"

#include "bench_matvec_fill_header.h"

int bench_matvec_fill_show_help(void){
  return 0;
}

const char* bench_matvec_fill_name(int argc, char** argv){
  (void)argc;
  (void)argv;

  static char name_buffer[MAX_LINE];
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [sgemv]");
  return name_buffer;
}

int bench_matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillX, double RealScaleX, double ImagScaleX, int incX, double RealBeta, double ImagBeta, int FillY, double RealScaleY, double ImagScaleY, int incY, int trials){
  int rc = 0;
  int i = 0;

  util_random_seed();
  int NX;
  int NY;
  switch(TransA){
    case 'n':
    case 'N':
      NX = N;
      NY = M;
    break;
    default:
      NX = M;
      NY = N;
    break;
  }

  float *A  = util_smat_alloc(Order, M, N, lda);
  float *X  = util_svec_alloc(NX, incX);
  float *Y  = util_svec_alloc(NY, incY);
  float alpha = RealAlpha;
  float beta = RealBeta;

  util_smat_fill(Order, 'n', M, N, A, lda, FillA, RealScaleA, ImagScaleA);
  util_svec_fill(NX, X, incX, FillX, RealScaleX, ImagScaleX);
  util_svec_fill(NY, Y, incY, FillY, RealScaleY, ImagScaleY);
  float *res  = (float*)malloc(NY * incY * sizeof(float));

  for(i = 0; i < trials; i++){
    memcpy(res, Y, NY * incY * sizeof(float));
    time_tic();
    CALL_SGEMV(Order, TransA, M, N, alpha, A, lda, X, incX, beta, res, incY);
    time_toc();
  }

  double dN = (double)N;
  double dM = (double)M;
  double dNY = (double)NY;
  metric_load_double("time", time_read());
  metric_load_double("trials", (double)(trials));
  metric_load_double("input", dN * dM + dN + dM);
  metric_load_double("output", dNY);
  metric_load_double("d_fma", dN * dM);
  metric_dump();

  free(X);
  free(Y);
  free(res);
  return rc;
}
