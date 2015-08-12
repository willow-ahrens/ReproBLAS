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

static void bench_rzgemv_options_initialize(void){
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
  bench_rzgemv_options_initialize();

  opt_show_option(fold);
  return 0;
}

const char* bench_matvec_fill_name(int argc, char** argv){
  (void)argc;
  (void)argv;
  bench_rzgemv_options_initialize();

  opt_eval_option(argc, argv, &fold);

  static char name_buffer[MAX_LINE];
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [rzgemv] (fold = %d)", fold._int.value);
  return name_buffer;
}

int bench_matvec_fill_test(int argc, char** argv, char Order, char TransA, int M, int N, double RealAlpha, double ImagAlpha, int FillA, double RealScaleA, double ImagScaleA, int lda, int FillX, double RealScaleX, double ImagScaleX, int incX, double RealBeta, double ImagBeta, int FillY, double RealScaleY, double ImagScaleY, int incY, int trials){
  int rc = 0;
  int i = 0;
  int j = 0;

  bench_rzgemv_options_initialize();

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

  double complex *A  = util_zmat_alloc(Order, M, N, lda);
  double complex *X  = util_zvec_alloc(NX, incX);
  double complex *Y  = util_zvec_alloc(NY, incY);
  double complex alpha = RealAlpha + I * ImagBeta;
  double complex beta = RealBeta + I * ImagBeta;
  double complex betaY;
  double_complex_indexed *YI = (double_complex_indexed*)malloc(NY * idxd_zisize(fold._int.value));

  util_zmat_fill(Order, 'n', M, N, A, lda, FillA, RealScaleA, ImagScaleA);
  util_zvec_fill(NX, X, incX, FillX, RealScaleX, ImagScaleX);
  util_zvec_fill(NY, Y, incY, FillY, RealScaleY, ImagScaleY);
  double complex *res  = (double complex*)malloc(NY * incY * sizeof(double complex));

  if(fold._int.value == DIDEFAULTFOLD){
    for(i = 0; i < trials; i++){
      memcpy(res, Y, NY * incY * sizeof(double complex));
      time_tic();
      rzgemv(Order, TransA, M, N, &alpha, A, lda, X, incX, &beta, res, incY);
      time_toc();
    }
  }else{
    for(i = 0; i < trials; i++){
      memcpy(res, Y, NY * incY * sizeof(double complex));
      time_tic();
      if(beta == 1.0){
        for(j = 0; j < NY; j++){
          idxd_zizconv(fold._int.value, res + j * incY, YI + j * idxd_zisize(fold._int.value));
        }
      }else{
        for(j = 0; j < NY; j++){
          betaY = res[j * incY] * beta;
          idxd_zizconv(fold._int.value, &betaY, YI + j * idxd_zisize(fold._int.value));
        }
      }
      zizgemv(fold._int.value, Order, TransA, M, N, &alpha, A, lda, X, incX, YI, 1);
      for(j = 0; j < NY; j++){
        idxd_zziconv_sub(fold._int.value, YI + j * idxd_zisize(fold._int.value), res + j * incY);
      }
      time_toc();
    }
  }

  metric_load_double("time", time_read());
  metric_load_double("trials", (double)(trials));
  metric_load_double("input", (double)(N * M + N + M));
  metric_load_double("output", (double)NY);
  metric_load_double("d_mul", (double)(4 * N * M));
  metric_load_double("d_add", (double)(3 * fold._int.value - 2) * 4 * N * M);
  metric_load_double("d_orb", (double)(4 * fold._int.value * N * M));
  metric_dump();

  free(X);
  free(Y);
  free(YI);
  free(res);
  return rc;
}
