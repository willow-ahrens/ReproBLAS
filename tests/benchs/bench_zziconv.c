#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <indexedBLAS.h>

#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_metric.h"

#include "../../config.h"

#include "bench_vecvec_fill_header.h"

static opt_option fold;

static void bench_zziconv_options_initialize(void){
  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 0;
  fold._int.max               = DIMAXFOLD;
  fold._int.value             = DIDEFAULTFOLD;
}

int bench_vecvec_fill_show_help(void){
  bench_zziconv_options_initialize();

  opt_show_option(fold);

  return 0;
}

const char* bench_vecvec_fill_name(int argc, char** argv){
  (void)argc;
  (void)argv;
  static char name_buffer[MAX_LINE];

  bench_zziconv_options_initialize();
  opt_eval_option(argc, argv, &fold);

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [zziconv] (fold = %d)", fold._int.value);
  return name_buffer;
}

int bench_vecvec_fill_test(int argc, char** argv, int N, int FillX, double RealScaleX, double ImagScaleX, int incX, int FillY, double RealScaleY, double ImagScaleY, int incY, int trials){
  (void)argc;
  (void)argv;
  (void)FillY;
  (void)RealScaleY;
  (void)ImagScaleY;
  (void)incY;
  int rc = 0;
  int i;
  int j;
  double complex res = 0.0;
  double_complex_indexed *ires;

  bench_zziconv_options_initialize();
  opt_eval_option(argc, argv, &fold);

  util_random_seed();

  double complex *X = util_zvec_alloc(N, incX);

  //fill x
  util_zvec_fill(N, X, incX, FillX, RealScaleX, ImagScaleX);

  if(fold._int.value == 0){
    ires = zialloc(DIMAXFOLD);
    zisetzero(DIMAXFOLD, ires);
    zizsum(DIMAXFOLD, N, X, incX, ires);
    time_tic();
    for(j = 1; j <= DIMAXFOLD; j++){
      for(i = 0; i < trials; i++){
        zziconv_sub(j, ires, &res);
      }
    }
    time_toc();
    free(ires);
  }else{
    ires = zialloc(fold._int.value);
    zisetzero(fold._int.value, ires);
    zizsum(fold._int.value, N, X, incX, ires);
    time_tic();
    for(i = 0; i < trials; i++){
      zziconv_sub(fold._int.value, ires, &res);
    }
    time_toc();
    free(ires);
  }

  metric_load_double("time", time_read());
  metric_load_double("res_real", crealf(res));
  metric_load_double("res_imag", cimagf(res));
  metric_load_double("trials", (double)trials);
  if(fold._int.value == 0){
    metric_load_double("input", (double)DIMAXFOLD);
  }else{
    metric_load_double("input", (double)1);
  }
  metric_load_double("output", (double)1);
  metric_dump();

  free(X);
  return rc;
}
