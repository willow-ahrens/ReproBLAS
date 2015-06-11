#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <indexedBLAS.h>
#include <reproBLAS.h>

#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_metric.h"

#include "../../config.h"

#include "bench_vecvec_fill_header.h"

static opt_option fold;

static void bench_rzsum_options_initialize(void){
  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 0;
  fold._int.max               = MAX_FOLD;
  fold._int.value             = DEFAULT_FOLD;
}

int bench_vecvec_fill_show_help(void){
  bench_rzsum_options_initialize();

  opt_show_option(fold);

  return 0;
}

const char* bench_vecvec_fill_name(int argc, char** argv){
  (void)argc;
  (void)argv;
  static char name_buffer[MAX_LINE];

  bench_rzsum_options_initialize();
  opt_eval_option(argc, argv, &fold);

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [rzsum] (fold = %d)", fold._int.value);
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

  bench_rzsum_options_initialize();
  opt_eval_option(argc, argv, &fold);

  util_random_seed();

  double complex *X = util_zvec_alloc(N, incX);

  //fill X
  util_zvec_fill(N, X, incX, FillX, RealScaleX, ImagScaleX);

  if(fold._int.value == DEFAULT_FOLD){
    time_tic();
    for(i = 0; i < trials; i++){
      rzsum_sub(N, X, incX, &res);
    }
    time_toc();
  }else if(fold._int.value == 0){
    time_tic();
    for(j = 1; j <= MAX_FOLD; j++){
      ires = zialloc(j);
      zisetzero(j, ires);
      for(i = 0; i < trials; i++){
        zizsum(j, N, X, incX, ires);
      }
      zziconv_sub(j, ires, &res);
      free(ires);
    }
    time_toc();
  }else{
    time_tic();
    ires = zialloc(fold._int.value);
    zisetzero(fold._int.value, ires);
    for(i = 0; i < trials; i++){
      zizsum(fold._int.value, N, X, incX, ires);
    }
    zziconv_sub(fold._int.value, ires, &res);
    free(ires);
    time_toc();
  }

  metric_load_double("time", time_read());
  metric_load_double("res_real", creal(res));
  metric_load_double("res_imag", cimag(res));
  metric_load_double("trials", (double)trials);
  metric_load_double("input", (double)N);
  metric_load_double("output", (double)1);
  if(fold._int.value == 0){
    metric_load_double("d_add", (double)(3 * MAX_FOLD * (MAX_FOLD - 1) * 0.5 - 2 * MAX_FOLD) * 2 * N);
    metric_load_double("d_orb", (double)MAX_FOLD * (MAX_FOLD - 1) * 0.5 * 2 * N);
  }else{
    metric_load_double("d_add", (double)(3 * fold._int.value - 2) * 2 * N);
    metric_load_double("d_orb", (double)fold._int.value * 2 * N);
  }
  metric_dump();

  free(X);
  return rc;
}
