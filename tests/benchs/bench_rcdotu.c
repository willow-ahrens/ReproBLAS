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

static void bench_rcdotu_options_initialize(void){
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
  bench_rcdotu_options_initialize();

  opt_show_option(fold);

  return 0;
}

const char* bench_vecvec_fill_name(int argc, char** argv){
  (void)argc;
  (void)argv;
  static char name_buffer[MAX_LINE];

  bench_rcdotu_options_initialize();
  opt_eval_option(argc, argv, &fold);

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [rcdotu] (fold = %d)", fold._int.value);
  return name_buffer;
}

int bench_vecvec_fill_test(int argc, char** argv, int N, int FillX, double RealScaleX, double ImagScaleX, int incX, int FillY, double RealScaleY, double ImagScaleY, int incY, int trials){
  (void)argc;
  (void)argv;
  int rc = 0;
  int i;
  int j;
  float complex res = 0.0;
  float_complex_indexed *ires;

  bench_rcdotu_options_initialize();
  opt_eval_option(argc, argv, &fold);

  util_random_seed();

  float complex *X = util_cvec_alloc(N, incX);
  float complex *Y = util_cvec_alloc(N, incY);

  //fill X and Y
  util_cvec_fill(N, X, incX, FillX, RealScaleX, ImagScaleX);
  util_cvec_fill(N, Y, incY, FillY, RealScaleY, ImagScaleY);

  if(fold._int.value == DEFAULT_FOLD){
    time_tic();
    for(i = 0; i < trials; i++){
      rcdotu_sub(N, X, incX, Y, incY, &res);
    }
    time_toc();
  }else if(fold._int.value == 0){
    time_tic();
    for(j = 1; j <= MAX_FOLD; j++){
      ires = cialloc(j);
      cisetzero(j, ires);
      for(i = 0; i < trials; i++){
        cicdotu(j, N, X, incX, Y, incY, ires);
      }
      cciconv_sub(j, ires, &res);
      free(ires);
    }
    time_toc();
  }else{
    time_tic();
    ires = cialloc(fold._int.value);
    cisetzero(fold._int.value, ires);
    for(i = 0; i < trials; i++){
      cicdotu(fold._int.value, N, X, incX, Y, incY, ires);
    }
    cciconv_sub(fold._int.value, ires, &res);
    free(ires);
    time_toc();
  }

  metric_load_double("time", time_read());
  metric_load_float("res_real", crealf(res));
  metric_load_float("res_imag", cimagf(res));
  metric_load_double("trials", (double)trials);
  metric_load_double("input", (double)2 * N);
  metric_load_double("output", (double)1);
  metric_load_double("s_mul", (double)4 * N);
  if(fold._int.value == 0){
    metric_load_double("s_add", (double)(3 * MAX_FOLD * (MAX_FOLD - 1) * 0.5 - 2 * MAX_FOLD) * 4 * N);
    metric_load_double("s_orb", (double)MAX_FOLD * (MAX_FOLD - 1) * 0.5 * 4 * N);
  }else{
    metric_load_double("s_add", (double)(3 * fold._int.value - 2) * 4 * N);
    metric_load_double("s_orb", (double)fold._int.value * 4 * N);
  }
  metric_dump();

  free(X);
  free(Y);
  return rc;
}
