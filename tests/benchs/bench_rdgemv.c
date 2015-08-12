#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <indexedBLAS.h>
#include <reproBLAS.h>

#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_metric.h"

#include "../../config.h"

#include "bench_matvec_fill_header.h"

static opt_option fold;

static void bench_rdgemv_options_initialize(void){
  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 0;
  fold._int.max               = DIMAXFOLD;
  fold._int.value             = DIDEFAULTFOLD;
}

int bench_matvec_fill_show_help(void){
  bench_rdgemv_options_initialize();

  opt_show_option(fold);
  return 0;
}

const char* bench_matvec_fill_name(int argc, char** argv){
  (void)argc;
  (void)argv;
  bench_rdgemv_options_initialize();

  opt_eval_option(argc, argv, &fold);

  static char name_buffer[MAX_LINE];
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [rdgemv] (fold = %d)", fold._int.value);
  return name_buffer;
}

int bench_matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillX, double RealScaleX, double ImagScaleX, int incX, double RealBeta, double ImagBeta, int FillY, double RealScaleY, double ImagScaleY, int incY, int trials){
  int rc = 0;
  int t;
  int i;
  int j;

  bench_rdgemv_options_initialize();

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

  double *A  = util_dmat_alloc(Order, M, N, lda);
  double *X  = util_dvec_alloc(NX, incX);
  double *Y  = util_dvec_alloc(NY, incY);
  double alpha = RealAlpha;
  double beta = RealBeta;
  double_indexed *YI = (double_indexed*)malloc(NY * disize(fold._int.value));

  util_dmat_fill(Order, 'n', M, N, A, lda, FillA, RealScaleA, ImagScaleA);
  util_dvec_fill(NX, X, incX, FillX, RealScaleX, ImagScaleX);
  util_dvec_fill(NY, Y, incY, FillY, RealScaleY, ImagScaleY);
  double *res  = (double*)malloc(NY * incY * sizeof(double));

  if(fold._int.value == DIDEFAULTFOLD){
    for(t = 0; t < trials; t++){
      memcpy(res, Y, NY * incY * sizeof(double));
      time_tic();
      rdgemv(Order, TransA, M, N, alpha, A, lda, X, incX, beta, res, incY);
      time_toc();
    }
  }else{
    for(t = 0; t < trials; t++){
      memcpy(res, Y, NY * incY * sizeof(double));
      time_tic();
      if(beta == 1.0){
        for(j = 0; j < NY; j++){
          didconv(fold._int.value, res[j * incY], YI + j * disize(fold._int.value));
        }
      }else{
        for(j = 0; j < NY; j++){
          didconv(fold._int.value, res[j * incY] * beta, YI + j * disize(fold._int.value));
        }
      }
      didgemv(fold._int.value, Order, TransA, M, N, alpha, A, lda, X, incX, YI, 1);
      for(j = 0; j < NY; j++){
        res[j * incY] = ddiconv(fold._int.value, YI + j * disize(fold._int.value));
      }
      time_toc();
    }
  }

  metric_load_double("time", time_read());
  metric_load_double("trials", (double)(trials));
  metric_load_double("input", (double)(N * M + N + M));
  metric_load_double("output", (double)NY);
  metric_load_double("d_mul", (double)(N * M));
  metric_load_double("d_add", (double)((3 * fold._int.value - 2) * N * M));
  metric_load_double("d_orb", (double)(fold._int.value * N * M));
  metric_dump();

  free(X);
  free(Y);
  free(YI);
  free(res);
  return rc;
}
