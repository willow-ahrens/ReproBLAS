#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <idxd.h>
#include <reproBLAS.h>

#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_metric.h"

#include "../../config.h"

#include "bench_matvec_fill_header.h"

static opt_option fold;

static void bench_rsgemv_options_initialize(void){
  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 2;
  fold._int.max               = idxd_DIMAXFOLD;
  fold._int.value             = DIDEFAULTFOLD;
}

int bench_matvec_fill_show_help(void){
  bench_rsgemv_options_initialize();

  opt_show_option(fold);
  return 0;
}

const char* bench_matvec_fill_name(int argc, char** argv){
  (void)argc;
  (void)argv;
  bench_rsgemv_options_initialize();

  opt_eval_option(argc, argv, &fold);

  static char name_buffer[MAX_LINE];
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [rsgemv] (fold = %d)", fold._int.value);
  return name_buffer;
}

int bench_matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillX, double RealScaleX, double ImagScaleX, int incX, double RealBeta, double ImagBeta, int FillY, double RealScaleY, double ImagScaleY, int incY, int trials){
  int rc = 0;
  int i;

  bench_rsgemv_options_initialize();

  opt_eval_option(argc, argv, &fold);

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
  float *res = util_svec_alloc(NY, incY);
  float alpha = RealAlpha;
  float beta = RealBeta;

  util_smat_fill(Order, 'n', M, N, A, lda, FillA, RealScaleA, ImagScaleA);
  util_svec_fill(NX, X, incX, FillX, RealScaleX, ImagScaleX);
  util_svec_fill(NY, Y, incY, FillY, RealScaleY, ImagScaleY);

  for(i = 0; i < trials; i++){
    memcpy(res, Y, NY * incY * sizeof(float));
    time_tic();
    reproBLAS_rsgemv(fold._int.value, Order, TransA, M, N, alpha, A, lda, X, incX, beta, res, incY);
    time_toc();
  }

  double dM = (double)M;
  double dN = (double)N;
  double dNY = (double)NY;
  metric_load_double("time", time_read());
  metric_load_double("trials", (double)trials);
  metric_load_double("input", dN * dM + dN + dM);
  metric_load_double("output", dNY);
  metric_load_double("normalizer", dN * dM);
  metric_load_double("s_mul", dN * dM);
  metric_load_double("s_add", (3 * fold._int.value - 2) * dN * dM);
  metric_load_double("s_orb", fold._int.value * dN * dM);
  metric_dump();

  free(A);
  free(X);
  free(Y);
  free(res);
  return rc;
}
