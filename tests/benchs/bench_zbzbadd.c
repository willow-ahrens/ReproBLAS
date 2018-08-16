#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <binnedBLAS.h>
#include <reproBLAS.h>

#include "../common/test_opt.h"
#include "../common/test_time.h"
#include "../common/test_metric.h"

#include "../../config.h"

#include "bench_vecvec_fill_header.h"

static opt_option fold;
static opt_option preN;

static void bench_zbzbadd_options_initialize(void){
  fold._int.header.type       = opt_int;
  fold._int.header.short_name = 'k';
  fold._int.header.long_name  = "fold";
  fold._int.header.help       = "fold";
  fold._int.required          = 0;
  fold._int.min               = 2;
  fold._int.max               = binned_DBMAXFOLD;
  fold._int.value             = DIDEFAULTFOLD;

  preN._int.header.type       = opt_int;
  preN._int.header.short_name = 'k';
  preN._int.header.long_name  = "preN";
  preN._int.header.help       = "zbzbadd preN before zbzbadd";
  preN._int.required          = 0;
  preN._int.min               = 1;
  preN._int.max               = INT_MAX;
  preN._int.value             = 1024;
}

int bench_vecvec_fill_show_help(void){
  bench_zbzbadd_options_initialize();

  opt_show_option(fold);

  return 0;
}

const char* bench_vecvec_fill_name(int argc, char** argv){
  (void)argc;
  (void)argv;
  static char name_buffer[MAX_LINE];

  bench_zbzbadd_options_initialize();
  opt_eval_option(argc, argv, &fold);

  snprintf(name_buffer, MAX_LINE * sizeof(char), "Benchmark [zbzbadd] (fold = %d)", fold._int.value);
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
  int k;
  double complex res = 0.0;
  double_complex_binned *ires;

  bench_zbzbadd_options_initialize();
  opt_eval_option(argc, argv, &fold);

  util_random_seed();

  double complex *preX = util_zvec_alloc(N * preN._int.value, incX);
  double_complex_binned *X;

  //fill x
  util_zvec_fill(N * preN._int.value, preX, incX, FillX, RealScaleX, ImagScaleX);

  X = (double_complex_binned*)util_zvec_alloc(N * binned_zbnum(fold._int.value), 1);
  memset(X, 0, N * binned_zbsize(fold._int.value));
  for(i = 0; i < N; i++){
    binnedBLAS_zbzsum(fold._int.value, preN._int.value, preX + i * preN._int.value * incX, incX, X + i * binned_zbnum(fold._int.value));
  }
  ires = binned_zballoc(fold._int.value);
  time_tic();
  for(i = 0; i < trials; i++){
    binned_zbsetzero(fold._int.value, ires);
    for(k = 0; k < N; k++){
      binned_zbzbadd(fold._int.value, X + k * binned_zbnum(fold._int.value), ires);
    }
    binned_zzbconv_sub(fold._int.value, ires, &res);
  }
  time_toc();
  free(ires);
  free(X);

  double dN = (double)N;
  metric_load_double("time", time_read());
  metric_load_double("res_real", creal(res));
  metric_load_double("res_imag", cimag(res));
  metric_load_double("trials", (double)trials);
  metric_load_double("input", dN);
  metric_load_double("output", 1.0);
  metric_load_double("normalizer", dN);
  metric_dump();

  free(preX);
  return rc;
}
